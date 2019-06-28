#!/usr/bin/env python

import copy
import json

def gen_qcfractal_geo_opt_json(input_json, qm_method, qm_basis, output_json):
    add_procedure_args = []
    options = {
        "keywords": None,
        "qc_spec": {
            "driver": "gradient",
            "method": qm_method,
            "basis": qm_basis,
            "program": "psi4"
        },
    }
    with open(input_json) as infile:
        molecule_data_list = json.load(infile)
    for mdata in molecule_data_list:
        initial_molecules = mdata['initial_molecules']
        cmiles_ids = mdata['cmiles_identifiers']
        for initial_molecule in initial_molecules:
            this_option = copy.deepcopy(options)
            this_option['cmiles_ids'] = cmiles_ids
            add_procedure_args.append(("optimization", "geometric", this_option, initial_molecule))
    with open(output_json, 'w') as outfile:
        json.dump(add_procedure_args, outfile, indent=2)

def main():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("in_json", help="Input json file that contains a list of {initial_molecules,cmiles_identifiers}")
    parser.add_argument("-m", "--qm_method", default='b3lyp-d3bj', help="QM method to use for optimization")
    parser.add_argument("-b", "--qm_basis", default='dzvp', help="QM basis to use for optimization")
    parser.add_argument("-o", "--output", default='geometry_opt_submit.json', help='Output JSON file name')
    args = parser.parse_args()

    gen_qcfractal_geo_opt_json(args.in_json, args.qm_method, args.qm_basis, args.output)

if __name__ == "__main__":
    main()