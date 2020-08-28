#!/usr/bin/env python

import os
import sys
import copy
import yaml
import json
import numpy as np

from forcebalance.molecule import Molecule, Elements, bohr2ang
import qcfractal.interface as ptl

from process_molecules import read_sdf_to_fb_mol, generate_torsion_index
from find_dihedrals import DihedralSelector
from torsion_submitter import TorsionSubmitter


def get_group1_1d_dihedrals(filename):
    m = read_sdf_to_fb_mol(filename)
    dihedral_selector = DihedralSelector(m)
    dihedral_filters = ['heavy_atoms', 'no_ring', 'unique_center_bond']
    dihedral_list = dihedral_selector.find_dihedrals(dihedral_filters=dihedral_filters)
    return dihedral_list

def submit_group1_1d(filenames, scan_conf_file, client_conf_file, to_json):
    submitter = TorsionSubmitter(scan_conf_file=scan_conf_file, client_conf_file=client_conf_file)
    for f in filenames:
        print(f"\n*** Submitting 1-D torsion scans for {f} ***")
        dihedral_list = get_group1_1d_dihedrals(f)
        submitter.submit_1d(f, dihedral_list)
    submitter.write_checkpoint()

def prepare_group1_1d_json(filenames, scan_conf_file):
    submitter = TorsionSubmitter(scan_conf_file=scan_conf_file)
    for f in filenames:
        print(f"\n*** Preparing 1-D torsion scans as JSON for {f} ***")
        dihedral_list = get_group1_1d_dihedrals(f)
        submitter.prepare_1d_json(f, dihedral_list)
    submitter.write_checkpoint()
    submitter.write_submitted_json("submit_torsion_options.json")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='*', help='Input sdf file for a single molecule')
    parser.add_argument("-s", "--scan_config", help='File containing configuration of QM scans')
    parser.add_argument("-c", "--client_config", help='File containing configuration of QCFractal Client')
    parser.add_argument("-j", "--save_json", action="store_true", default=False, help='If specified, will not submit but save all submit job in a json.')
    args = parser.parse_args()

    print(' '.join(sys.argv))
    if args.save_json:
        prepare_group1_1d_json(args.infiles, args.scan_config)
    else:
        submit_group1_1d(args.infiles, args.scan_config, args.client_config)

if __name__ == '__main__':
    main()
