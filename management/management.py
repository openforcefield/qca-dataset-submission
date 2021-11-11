from time import sleep

from collections import defaultdict, Counter
from pprint import pformat

import qcportal as ptl


## failures and incompletes

### basic results

def get_results(dataset, method, basis, program):
    return dataset.get_records(method=method, basis=basis, program=program)


def get_unfinished_results(dataset, method, basis, program):
    res = get_results(dataset, method, basis, program)
    res = [r for r in res if r.status != 'COMPLETE']

    return res

### optimizations

def get_optimizations(dataset, spec, dropna=False):
    ds = dataset
    ds.status(spec)

    if dropna:
        df = ds.df.dropna()
    else:
        df = ds.df

    res = df[spec].tolist()
    return res


def get_unfinished_optimizations(dataset, spec, client, dropna=False):
    res = get_optimizations(dataset, spec, client, dropna=dropna)
    res = [opt for opt in res if opt.status != 'COMPLETE']

    return res


def _query_procedures(ids, client):
    res = []
    ids = list(ids)
    for i in range(0,len(ids),1000):
        ids_i = ids[i:i+1000]
        res_i = client.query_procedures(ids_i)
        res.extend(res_i)
        
    return res

### torsiondrives

def get_torsiondrives(dataset, spec):

    ds = dataset

    while True:
        try:
            ds.status(spec)
        except:
            pass
        else:
            break

    return ds.df[spec].tolist()


def get_torsiondrive_optimizations(
        dataset, spec, client, noncomplete=False):

    ds = dataset

    while True:
        try:
            ds.status(spec)
        except:
            pass
        else:
            break

    optimizations = defaultdict(set)
    for tdr in ds.df[spec]:
        for val in tdr.optimization_history.values():
            optimizations[tdr.id].update(set(val))

    res_opt = {key: _query_procedures(value, client) for key, value in optimizations.items()}
    if noncomplete:
        res_opt = {key: [opt for opt in value if opt.status != 'COMPLETE']
                   for key, value in res_opt.items()}

    return res_opt


def get_unfinished_torsiondrive_optimizations(
        dataset, spec, client, noncomplete=False):

    res_opt = get_torsiondrive_optimizations(dataset, spec, client, noncomplete=noncomplete)

    ds = dataset
    for tdr in ds.df[spec]:
        if tdr.status == 'COMPLETE':
            res_opt.pop(tdr.id, None)

    return res_opt


def merge(datadict):
    res_s = set()
    res = list()
    for val in datadict.values():
        new_ids = set(i.id for i in val) - res_s
        res.extend([i for i in val if i.id in new_ids])
    return res


## error messages

def count_unique_result_error_messages(
        results, client, full=False, pretty_print=False, tolerate_missing=False):
    errors = defaultdict(set)

    for res in results:
        if not isinstance(res, dict):
            r = res.dict()
        else:
            r = res

        if r['status'] != 'ERROR':
            continue

        try:
            kv = client.query_kvstore([r['id']])[r['id']]
            err_content = kv.get_json()
        except:
            err_content = None

        if tolerate_missing:
            if err_content is None:
                errors[None].add(r['id'])
                continue

        if full:
            errors[err_content['error_message']].add(r['id'])
        else:
            errors[err_content['error_message'].split('\n')[-2]].add(r['id'])

    errors = dict(errors)

    # Convert None key to a string to avoid issues with sorting below
    if None in errors:
        errors["None"] = errors[None]
        errors.pop(None)

    content = ""
    if pretty_print:
        for count, key, value in sorted([(len(value), key, value) for key, value in errors.items()], reverse=True):
            content += '-------------------------------------\n'
            content += f"count : {count}\n"
            content += '\n'
            content += f'{key}\n'
            content += '\n'
            content += 'ids : \n'
            content += f'{pformat(value, width=80, compact=True)}\n'
            content += '-------------------------------------\n'
        return content
    else:
        return errors


def get_unique_optimization_error_messages(optimizations, full=False):
    if full:
        return set(opt.get_error().error_message
                   for opt in optimizations if opt.status == 'ERROR')
    else:
        return set(opt.get_error().error_message.split('\n')[-2]
                   for opt in optimizations if opt.status == 'ERROR')


count_unique_optimization_error_messages = count_unique_result_error_messages


## restarts

def restart_results(results, client):
    erred = []
    for res in results:
        if not isinstance(res, dict):
            if res.status == 'ERROR':
                erred.append(res.id)
        else:
            if res['status'] == 'ERROR':
                erred.append(res['id'])

    client.modify_tasks(operation='restart', base_result=erred)
    print(f"Restarted ERRORed results:\n====\n{erred}\n====\n")


def regenerate_results(results, client):
    regenerate = []
    for res in results:
        if not isinstance(res, dict):
            if res.status == 'INCOMPLETE':
                regenerate.append(res.id)
        else:
            if res['status'] == 'INCOMPLETE':
                regenerate.append(res['id'])

    for rid in regenerate:
        print(f"Regnerated INCOMPLETE result `{rid}`")
        client.modify_tasks(operation='regenerate', base_result=rid)


restart_optimizations = restart_results


def regenerate_optimizations(optimizations, client):

    regenerate = []
    for res in optimizations:
        if not isinstance(res, dict):
            if res.status == 'INCOMPLETE' and (res.final_molecule is not None):
                regenerate.append(res.id)
        else:
            if res['status'] == 'INCOMPLETE' and (res['final_molecule'] is not None):
                regenerate.append(res['id'])

    for optid in regenerate:
        print(f"Regnerated INCOMPLETE optimization `{optid}`")
        client.modify_tasks(operation='regenerate', base_result=optid)


def restart_torsiondrives(torsiondrives, client):
    for tdr in torsiondrives:
        if tdr.status == 'ERROR':
            print(f"Restarted ERRORed torsiondrive `{tdr.id}`")
            client.modify_services('restart', procedure_id=tdr.id)


## modify tasks

def retag_results(results, client, compute_tag):
    retag = []
    for res in results:
        if not isinstance(res, dict):
            if res.status != 'COMPLETE':
                retag.append(res.id)
        else:
            if res['status'] != 'COMPLETE':
                retag.append(res['id'])

    client.modify_tasks(operation='modify', base_result=retag, new_tag=compute_tag)
    print(f"Retagged results with `{compute_tag}':\n====\n{retag}\n====\n")


retag_optimizations = retag_results


def reprioritize_results(results, client, priority):
    """Priority can be one of "high", "normal", "low".

    """
    from qcportal.models.task_models import PriorityEnum
    priority_map = {"high": PriorityEnum.HIGH,
                    2: PriorityEnum.HIGH,
                    "normal": PriorityEnum.NORMAL,
                    1: PriorityEnum.NORMAL,
                    "low": PriorityEnum.LOW,
                    0: PriorityEnum.LOW,
                    }

    reprioritize = []
    for res in results:
        if not isinstance(res, dict):
            if res.status != 'COMPLETE':
                reprioritize.append(res.id)
        else:
            if res['status'] != 'COMPLETE':
                reprioritize.append(res['id'])

    client.modify_tasks(operation='modify', base_result=reprioritize, new_priority=priority_map[priority])
    print(f"Reprioritized results as `{priority}':\n====\n{reprioritize}\n====\n")


reprioritize_optimizations = reprioritize_results
