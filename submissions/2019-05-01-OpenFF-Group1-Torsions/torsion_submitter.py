import os
import sys
import copy
import yaml
import json
import numpy as np

from forcebalance.molecule import Molecule, Elements, bohr2ang
import qcfractal.interface as ptl

from find_dihedrals import DihedralSelector
from process_molecules import read_sdf_to_fb_mol, generate_torsion_index


class TorsionSubmitter:
    def __init__(self, scan_conf_file=None, client_conf_file=None):
        # load scan config from scan_conf_file
        self.scan_conf = self.load_scan_conf(scan_conf_file)
        # create a client from config file
        self.client = ptl.FractalClient.from_file(client_conf_file) if client_conf_file is not None else None
        # load checkpoint file
        self.checkpoint_filename = "torsion_submit_checkpoint.json"
        self.load_checkpoint(self.checkpoint_filename)
        # sync scan config with checkpoint
        self.sync_scan_conf()
        # hold all json option to submit
        self.json_submit_options = []

    def load_checkpoint(self, filename):
        if os.path.isfile(filename):
            with open(filename) as infile:
                self.state = json.load(infile)
        else:
            print(f"checkpoint file {filename} not found, starting with empty state")
            self.state = {}

    def write_checkpoint(self):
        with open(self.checkpoint_filename, 'w') as outfile:
            json.dump(self.state, outfile, indent=2)

    def load_scan_conf(self, filename=None):
        """
        Get the scan configuration from a yaml file
        Parameters
        ----------
        filename: str or None
        The input scan config filename (in yaml format). If None, default conf will be used.
        Returns
        -------
        scan_conf: dict
        """
        required_scan_conf = {
            'qm_method': 'b3lyp',
            'qm_basis': '6-31G*',
            'grid_spacing': 15,
        }
        optional_scan_conf = {
            'energy_upper_limit': 0.05,
        }
        default_scan_conf = {**required_scan_conf, **optional_scan_conf}
        if filename is None:
            print("scan_configure file not provided, using default config")
            fn = "td_scan_configure.yaml"
            conf = copy.deepcopy(default_scan_conf)
            with open(fn, 'w') as outfile:
                yaml.dump(conf, outfile, default_flow_style=False)
            print(f"scan configure saved as {fn}")
            return conf
        with open(filename) as infile:
            conf = yaml.load(infile, Loader=yaml.SafeLoader)
        # convert keys to lower case
        conf = {k.lower():v for k,v in conf.items()}
        # check redundant and missing keys
        diff1 = required_scan_conf.keys() - conf.keys()
        if diff1:
            raise ValueError(f"Keys missing in scan_config file {filename}:\n {diff1}")
        diff2 = conf.keys() - default_scan_conf.keys()
        if diff2:
            print(f"Warning: Keys in scan_config file {filename} are ignored:\n {diff2}")
        return conf

    def sync_scan_conf(self):
        if 'scan_conf' in self.state:
            existing_conf = self.state['scan_conf']
            current_conf = self.scan_conf
            # make sure the current configuration matches the previous one
            assert json.dumps(existing_conf, sort_keys=True) == json.dumps(current_conf, sort_keys=True), \
                f'Error: previous scan conf {existing_conf} not consistent with current conf {current_conf}'
        else:
            self.state['scan_conf'] = self.scan_conf

    def fb_molecule_to_qc_molecule(self, fb_molecule):
        """ Convert an forcebalance.molecule.Molecule object to a qcportal.Molecule object"""
        moldata = {
            'symbols': fb_molecule.elem,
            'geometry': fb_molecule.xyzs[0] / bohr2ang,
            'molecular_charge': fb_molecule.Data.get('molecular_charge', 0),
            'molecular_multiplicity': fb_molecule.Data.get('mult', 1),
            'connectivity': [list(bond) + [bond_order] for bond, bond_order in zip(fb_molecule.bonds, fb_molecule.bond_orders)],
        }
        # append the cmiles_id as identifier
        #cmiles_id = fb_molecule.Data.get('cmiles_id')
        #if cmiles_id is not None:
            #for discard_key in ['provenance', 'unique_protomer_representation', 'unique_tautomer_representation']:
            #    cmiles_id.pop(discard_key, None)
            #moldata['identifiers'] = cmiles_id
        return ptl.Molecule.from_data(moldata, dtype="dict", units="angstrom")

    def qc_molecule_to_fb_molecule(self, qc_molecule):
        """ Convert an qcportal.Molecule object to a forcebalance.molecule.Molecule object"""
        m = Molecule()
        m.elem = [Elements[i] for i in qc_molecule.atomic_numbers]
        m.xyzs = [qc_molecule.geometry * bohr2ang]
        m.molecular_charge = qc_molecule.molecular_charge
        m.mult = qc_molecule.molecular_multiplicity
        return m

    def make_submit_options(self, initial_molecule, dihedrals, to_json=False):
        if isinstance(self.scan_conf['grid_spacing'], (list, tuple)):
            grid_spacing = copy.deepcopy(self.scan_conf['grid_spacing'])
        else:
            grid_spacing = [self.scan_conf['grid_spacing']]
        torsiondrive_options = {
            "initial_molecule": initial_molecule,
            "keywords": {
                "dihedrals": dihedrals,
                "grid_spacing": grid_spacing,
            },
            "optimization_spec": {
                "program": "geometric",
                "keywords": {
                    "coordsys": "tric",
                    "enforce": 0.1,
                    "reset": True,
                    "qccnv": True,
                    "epsilon": 0.0,
                }
            },
            "qc_spec": {
                "driver": "gradient",
                "method": self.scan_conf['qm_method'],
                "basis": self.scan_conf['qm_basis'],
                "keywords": None,
                "program": "psi4",
            },
        }
        if 'energy_upper_limit' in self.scan_conf:
            torsiondrive_options['keywords']["energy_upper_limit"] = self.scan_conf['energy_upper_limit']
        return torsiondrive_options

    def submit_1d(self, filename, dihedral_list):
        # read mol file and get cmiles id
        m = read_sdf_to_fb_mol(filename)
        qc_mol = self.fb_molecule_to_qc_molecule(m)
        mol_id = self.client.add_molecules([qc_mol])[0]
        # all options to be submitted
        all_job_options = []
        for dihedral in dihedral_list:
            torsiondrive_submit_option = self.make_submit_options(mol_id, [dihedral])
            all_job_options.append(torsiondrive_options)
        print(f"Submitting {len(all_job_options)} torsiondrive jobs")
        r = self.client.add_service(all_job_options)
        job_id_list = r.ids
        assert len(job_id_list) == len(dihedral_list)
        # save in checkpoint
        self.state[filename] = {'dihedrals': {}}
        for i, dihedral in enumerate(dihedral_list):
            dihedral_name = '-'.join(map(str, dihedral))
            dihedral_data = {'jobid': job_id_list[i], 'status': 'submitted'}
            self.state[filename]['dihedrals'][dihedral_name] = dihedral_data

    def prepare_1d_json(self, filename, dihedral_list):
        # read mol file and get cmiles id
        m = read_sdf_to_fb_mol(filename)
        qc_mol = self.fb_molecule_to_qc_molecule(m)
        mol_json = qc_mol.json_dict()
        cmiles_id = m.Data.get('cmiles_id', {})
        # all options to be submitted
        all_job_options = []
        canonical_torsion_labels = []
        for dihedral in dihedral_list:
            torsiondrive_submit_option = self.make_submit_options(mol_json, [dihedral])
            canonical_torsion_label = generate_torsion_index(m, dihedral)
            canonical_torsion_labels.append(canonical_torsion_label)
            # append the id information
            torsiondrive_submit_option["attributes"] = {
                'canonical_torsion_label': canonical_torsion_label,
                'cmiles_id': cmiles_id,
            }
            all_job_options.append(torsiondrive_submit_option)
        print(f"Saving {len(all_job_options)} torsiondrive options to json")
        self.json_submit_options.extend(all_job_options)
        # save in checkpoint
        self.state[filename] = {'dihedrals': {}}
        for i, dihedral in enumerate(dihedral_list):
            dihedral_name = '-'.join(map(str, dihedral))
            dihedral_data = {'status': 'saved_json', 'canonical_torsion_label': canonical_torsion_labels[i]}
            self.state[filename]['dihedrals'][dihedral_name] = dihedral_data

    def submit_2d(self, filename, dihedral_pairs_list):
        # read mol file and get cmiles id
        m = read_sdf_to_fb_mol(filename)
        qc_mol = self.fb_molecule_to_qc_molecule(m)
        mol_id = self.client.add_molecules([qc_mol])[0]
        # all options to be submitted
        all_job_options = []
        for dihedral_pair in dihedral_pairs_list:
            torsiondrive_submit_option = self.make_submit_options(mol_id, dihedral_pair)
            all_job_options.append(torsiondrive_options)
        print(f"Submitting {len(all_job_options)} torsiondrive jobs")
        r = self.client.add_service(all_job_options)
        job_id_list = r.ids
        assert len(job_id_list) == len(dihedral_list)
        # save in checkpoint
        self.state[filename] = {'dihedrals': {}}
        for i, dihedral_pair in enumerate(dihedral_pairs_list):
            d1, d2 = dihedral_pair
            dihedral_name = '-'.join(map(str, d1)) + '_' + '-'.join(map(str, d2))
            dihedral_data = {'jobid': job_id_list[i], 'status': 'submitted'}
            self.state[filename]['dihedrals'][dihedral_name] = dihedral_data

    def prepare_2d_json(self, filename, dihedral_pairs_list):
        # read mol file and get cmiles id
        m = read_sdf_to_fb_mol(filename)
        qc_mol = self.fb_molecule_to_qc_molecule(m)
        mol_json = qc_mol.json_dict()
        cmiles_id = m.Data.get('cmiles_id', {})
        # all options to be submitted
        all_job_options = []
        canonical_torsion_labels = []
        for dihedral1, dihedral2 in dihedral_pairs_list:
            torsiondrive_submit_option = self.make_submit_options(mol_json, [dihedral1, dihedral2])
            # append the id information
            canonical_torsion_label = generate_torsion_index(m, dihedral1) + ',' + generate_torsion_index(m, dihedral2)
            canonical_torsion_labels.append(canonical_torsion_label)
            torsiondrive_submit_option["attributes"] = {
                'canonical_torsion_label': canonical_torsion_label,
                'cmiles_id': cmiles_id,
            }
            all_job_options.append(torsiondrive_submit_option)
        print(f"Saving {len(all_job_options)} torsiondrive options to json")
        self.json_submit_options.extend(all_job_options)
        # save in checkpoint
        self.state[filename] = {'dihedrals': {}}
        for dihedral_pair in dihedral_pairs_list:
            d1, d2 = dihedral_pair
            dihedral_name = '-'.join(map(str, d1)) + '_' + '-'.join(map(str, d2))
            dihedral_data = {'status': 'saved_json', 'canonical_torsion_label': canonical_torsion_labels[i]}
            self.state[filename]['dihedrals'][dihedral_name] = dihedral_data

    def get_charge(self, filename):
        """ Get charge of molecule based on filename generated by OpenEye """
        name = os.path.splitext(filename)[0]
        if name[-1] == '+':
            charge = 1
        elif name[-1] == '-':
            charge = -1
        else:
            charge = 0
        return charge

    def write_submitted_json(self, filename):
        with open(filename, 'w') as outfile:
            json.dump(self.json_submit_options, outfile, indent=2)
            print(f'Total {len(self.json_submit_options)} torsiondrive options written to {filename}')