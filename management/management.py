from collections import defaultdict

import qcportal as ptl


def get_failed_optimizations(dataset, spec, client):
    ds = dataset

    while True:
        try:
            ds.status(spec)
        except:
            pass
        else:
            break

    ids = list(set(i.id for i in ds.df[spec]))
    
    res = []
    for i in range(0,len(ids),1000):
        ids_i = ids[i:i+1000]
        res_i = client.query_procedures(ids_i)
        res.extend(res_i)
        
    return res


def get_incomplete_torsiondrive_optimizations(dataset, spec, client, merged=False):
    ds = dataset

    while True:
        try:
            ds.status(spec)
        except:
            pass
        else:
            break

    ids = set(i.id for i in ds.df[spec])
    res = client.query_procedures(ids)
    
    optimizations = defaultdict(set)
    for tdr in ds.df.default:
        if tdr.status == 'COMPLETE':
            continue
            
        for val in tdr.optimization_history.values():
            optimizations[tdr.id].update(set(val))

    if merged:
        optimizations_i = set()
        for i in optimizations.values():
            optimizations_i.update(set(i))
        
        res_opt = client.query_procedures(optimizations_i)
    else:
        res_opt = {key: client.query_procedures(value) for key, value in optimizations.items()}

    return res_opt
