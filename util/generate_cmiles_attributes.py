#!/usr/bin/env python3
import json
import sys
import os
import io

import numpy as np

def dict_to_xyz_string(mol):
    bohr2ang = 0.529177
    out = ""
    out += "{:d}\n\n".format(len(mol['symbols']))
    for elem, xyz in zip(mol['symbols'], np.reshape(mol['geometry'], (-1,3))):
        out += "{:3s} {:8.6f} {:8.6f} {:8.6f}\n".format(elem, *[q*bohr2ang for q in xyz])
    return out


def qcmiles(jsmol, toolkit='rdkit'):

    import cmiles
    from openbabel import openbabel
    from openforcefield.topology.molecule import Molecule

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "sdf")
    obmol = openbabel.OBMol()

    xyz_str = dict_to_xyz_string(jsmol)

    obConversion.ReadString(obmol, xyz_str)
    sdf = obConversion.WriteString(obmol)

    with io.StringIO(sdf) as sdf_stream:
        qcmol = Molecule.from_file(sdf_stream, file_format='SDF').to_qcschema().dict()

    # would be nice if oFF could handle data as strings to avoid IO
    #with open('mol.sdf','w') as f:
    #    f.write(sdf)

    #qcmol = Molecule.from_file('mol.sdf').to_qcschema().dict()
    #os.remove("mol.sdf")

    # cmiles wants the schema with a flat xyz
    if len(qcmol['geometry'].shape) > 1:
        qcmol['geometry'] = np.reshape(qcmol['geometry'], (-1,))

    attribs = cmiles.generator.get_molecule_ids(qcmol, toolkit=toolkit)
    return attribs


def qcmiles_from_json_file(fnm, toolkit='rdkit'):

    with open(fnm, 'r') as f:
        js = json.load(f)

    return qcmiles_from_json_str(js, toolkit=toolkit)


def qcmiles_from_json_str(js, toolkit='rdkit'):

    js = js if isinstance(js, list) else [js]
    return [qcmiles(mol, toolkit=toolkit) for mol in js]


def main():

    import argparse

    parser = argparse.ArgumentParser(description="Generates CMILES data for use with QCArchive when connectivity is not known. Requires cmiles, the openff toolkit, OpenEye for protomers, RDkit for the other attributes, and openbabel")
    parser.add_argument("input_json", help="input file in json format. The json should be a list of QCSchema-style symbols and coordinates (in bohr). A dict(geometry=xyz, symbols=syms) is all that is needed for each molecule.")
    parser.add_argument("--use-rdkit", default=False, action="store_true", help="CMILES generator should use rdkit when possible")
    parser.add_argument("--use-openeye", default=False, action="store_true", help="CMILES generator should use only OpenEye")

    args = parser.parse_args().__dict__

    if args["use_openeye"] and args["use_rdkit"]:
        raise Exception("Cannot specifiy both use-rdkit and use-openeye")

    if not (args["use_rdkit"] or args["use_openeye"]):
        args["use_rdkit"] = True
    
    if args["use_rdkit"]:
        tk = 'rdkit'
    else:
        tk = 'openeye'

    print("using tk",tk)

    import pprint

    ret = qcmiles_from_json_file(args["input_json"], toolkit=tk)
    pp = pprint.PrettyPrinter(indent=4, width=80, compact=True)
    pp.pprint(ret)
    # for attr in ret:
    #     pp.pprint(attr)


if __name__ == "__main__":
    main()

