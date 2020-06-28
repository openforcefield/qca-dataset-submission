from collections import defaultdict, Counter

import qcportal as ptl


## failures and incompletes

def get_optimizations(dataset, spec, client, dropna=False):
    ds = dataset

    while True:
        try:
            ds.status(spec)
        except:
            pass
        else:
            break

    if dropna:
        df = ds.df.dropna()
    else:
        df = ds.df

    ids = list(set(i.id for i in df[spec]))
    
    res = _query_procedures(ids, client)    

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

def get_unique_optimization_error_messages(optimizations, full=False):
    if full:
        return set(opt.get_error().error_message
                   for opt in optimizations if opt.status == 'ERROR')
    else:
        return set(opt.get_error().error_message.split('\n')[-2]
                   for opt in optimizations if opt.status == 'ERROR')


def count_unique_optimization_error_messages(
        optimizations, full=False, pretty_print=False, tolerate_missing=False):
    errors = Counter()

    for opt in optimizations:
        if opt.status != 'ERROR':
            continue

        err_content = opt.get_error()

        if tolerate_missing:
            if err_content is None:
                errors += Counter({None: 1})
                continue

        if full:
            errors += Counter({err_content.error_message: 1})
        else:
            errors += Counter({err_content.error_message.split('\n')[-2]: 1})

    errors = dict(errors)

    if pretty_print:
        for key, value in errors.items():
            print(f"There are {value} instances of")
            print('-------------------------------------')
            print(key)
            print('-------------------------------------\n')
        return None
    else:
        return errors


## restarts

def restart_optimizations(optimizations, client):
    for opt in optimizations:
        if opt.status == 'ERROR':
            print(opt)
            client.modify_tasks(operation='restart', base_result=opt.id)


def restart_torsiondrives(torsiondrives, client):
    for tdr in torsiondrives:
        if tdr.status == 'ERROR':
            print(tdr)
            client.modify_services('restart', procedure_id=tdr.id)
