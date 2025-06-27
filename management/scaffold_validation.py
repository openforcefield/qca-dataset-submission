
import glob
import json
import copy
import warnings
from collections import defaultdict

from pprint import pprint

import numpy as np
import networkx as nx
import basis_set_exchange as bse
import psi4
from psi4.driver.qcdb import BasisSetNotFound

from openff.qcsubmit.constraints import Constraints
from openff.qcsubmit.exceptions import ConstraintError

psi4.core.be_quiet()

def check_basis_coverage(specifications, elements):
    """Check that all elements are supported by the basis set

    Parameters
    ----------
    specifications : dict
        Specifications as written from qcportal encode_to_json
    elements : list[str]
        List of element symbols

    """
    
    basis_report = defaultdict()
    for spec_name, spec in specifications.items():
        if "optimization_specification" in spec["specification"]: # TD
            spec = spec["specification"]["optimization_specification"]["qc_specification"]
        elif "qc_specification" in spec["specification"]: # Opt
            spec = spec["specification"]["qc_specification"]
        else: # SP or Opt
            spec = spec["specification"]
            
        if spec["program"].lower() == "torchani": # From QCSubmit _BaseDataset._get_missing_basis_coverage()
                # check ani1 first
                ani_coverage = {
                    "ani1x": {"C", "H", "N", "O"},
                    "ani1ccx": {"C", "H", "N", "O"},
                    "ani2x": {"C", "H", "N", "O", "S", "F", "Cl"},
                }
                covered_elements = ani_coverage[spec.method.lower()]
                # this is validated at the spec level so we should not get an error here
                difference = set(elements) - set(covered_elements)
        elif spec["program"].lower() == "psi4":
            basis = spec["basis"]
            if basis is None or basis.lower() == "none":
                warnings.warn(f"The spec {spec_name} has a basis of None, this will not be validated.")
                difference = set()
            elif "bse:" not in basis:
                # Use psi4.core.BasisSet to check if basis is available for each element
                difference = set()
                for elem in elements:
                    try:
                        psi4.core.BasisSet.build(psi4.geometry(f"{elem}"), "BASIS", basis);
                    except BasisSetNotFound:
                        difference.add(elem)
            else:
                basis_data = bse.get_basis(basis[4:], elements=elements, fmt='json')
                difference = set([elem for elem in elements if elem not in basis_data['elements']])
        elif spec["program"].lower() == "openmm":
            # smirnoff covered elements
            covered_elements = {"C", "H", "N", "O", "P", "S", "Cl", "Br", "F", "I"}
            difference = set(elements) - set(covered_elements)  
        elif spec["program"].lower() in ["xtb", "rdkit"]:
            difference = set()
        else:
            raise ValueError(f'Validation for {spec["program"].lower()} is not supported.')
                
        basis_report[spec_name] = difference
    
    return basis_report


def check_scf_keywords(keywords):
    """Check the SCF properties are valid oeprop keywords in psi4

    Returns
    -------
    bool
        True if all properties are valid, False otherwise.
    """
    keywords = copy.deepcopy(keywords)

    # Check Function Properties
    valid_props = {
        "DIPOLE_POLARIZABILITIES",
        "QUADRUPOLE_POLARIZABILITIES",
        "TRACELESS_QUADRUPOLE_POLARIZABILITIES",
        # Add more as needed from psi4 documentation
    }
    function_kwargs = keywords.pop("function_kwargs", {})
    if function_kwargs and "properties" in function_kwargs:
        for prop in function_kwargs["properties"]:
            try:
                if prop.upper() not in valid_props:
                    raise ValueError(f"PROP, {prop}, is not recognized.")
            except Exception as e:
                print(str(e))
                return False
        
    # Check SCF Properties
    # List of valid oeprop keywords in psi4 (as of psi4 v1.9+)
    valid_props = {
        "DIPOLE",
        "QUADRUPOLE",
        "ESP_AT_NUCLEI",
        "MULLIKEN_CHARGES",
        "LOWDIN_CHARGES",
        "LOWDIN_SPINS",
        "WIBERG_LOWDIN_INDICES",
        "NO_OCCUPATIONS",
        "MAYER_INDICES",
        "MBIS_CHARGES",
        "MBIS_VOLUME_RATIO",
        "MO_EXTENTS",
        # Add more as needed from psi4 documentation
    }
    scf_properties = keywords.pop("scf_properties")
    for prop in scf_properties:
        try:
            if prop.upper() not in valid_props and not prop.startswith("MULTIPOLE("):
                raise ValueError(f"OEPROP, {prop}, is not recognized.")
        except Exception as e:
            print(str(e))
            return False

        
    return True


def create_spec_report(spec, validated, extras):
    """Create specification report

    Parameters
    ----------
    spec : dict
        Specification dictionary
    validated : bool
        Whether the specification has been validated
    extras : dict
        Additional information to print

    Returns
    -------
    dict
        Report to print
    """

#    solvent = spec.get("implicit_solvent", None)
    spec_name = spec["name"]
    if "optimization_specification" in spec["specification"]: # TD
        spec = spec["specification"]["optimization_specification"]["qc_specification"]
    elif "qc_specification" in spec["specification"]: # Opt
        spec = spec["specification"]["qc_specification"]
    else: # SP or Opt
        spec = spec["specification"]

    if ("ddx" in spec["keywords"] or "pcm" in spec["keywords"]):
        # QCArchive passes implicit solvent keywords to Psi4 with keywords that start with the implicit 
        # solvent type, e.g., "specification": {"keyords": {"ddx_model": "pcm",...}...}. These keywords 
        # tend to start with the three letter names, 'pcm' or 'ddx', so we check for these prefixes.
        # https://psicode.org/psi4manual/master/autodoc_glossary_options_c.html#term-DDX_SOLVENT-DDX
        solvent = {key: val for key, val in spec["keywords"].items() if key[:3].lower() in ["pcm", "ddx"]}
    else:
        solvent = None

    data = {
        "**Specification Name**": spec_name,
        "**Method**": spec["method"],
        "**Basis**": spec["basis"],
        "**Wavefunction Protocol**": spec["protocols"].get("wavefunction", None),
        "**Implicit Solvent**": json.dumps(solvent) if solvent is not None else solvent,
        "**Keywords**": json.dumps({} if not spec.get("keywords", None) else spec["keywords"]),
        "**Validated**": validated
    }
    data.update(extras)
    return data

    
def get_constraint_checks(constraints, bonds):
    """Check that constraints are formatted correctly.

    Parameters
    ----------
    constraints : dict
        Constraint dictionary containing constraint types and information needed according to
        `geomeTRIC <https://geometric.readthedocs.io/en/latest/constraints.html#constraint-types>`_.
    bonds : list
        List of lists containing the indices of atoms connected to the atom represented by that index.

    Returns
    -------
    dict
        List of keys for which there are errors

    Raises
    ------
    ValueError
        Provided indices are not unique
    ConstraintError
        ``Constraint`` class from QCSubmit could not be initiated.
    """
    constraints_obj = Constraints()
    for con_type, constraint in constraints.items():
        try:
            if len(constraint["indices"]) != len(set(constraint["indices"])):
                raise ValueError("Constraints must have unique indices.")
            
            if con_type.lower() == "freeze":
                constraints_obj.add_freeze_constraint(
                    constraint_type=constraint["type"], indices=constraint["indices"], bonded=True
                )
            elif con_type.lower() == "set":
                constraints_obj.add_set_constraint(
                    constraint_type=constraint["type"], 
                    indices=constraint["indices"], 
                    bonded=True,
                    value=constraint["value"],
                )
            else:
                raise ConstraintError(
                    f"The constraint {constraint} is not available please chose from 'freeze' or 'set'."
                )
        except Exception as e:
            print(str(e))
            return {"constraint": True}

        pairs = [
            [constraint["indices"][i], constraint["indices"][i+1]] 
            for i in range(len(constraint["indices"]) - 1)
        ]
        try:
            for i, j in pairs:
                if j not in bonds[i]:
                    raise ValueError(f"Constraint pair ({i}, {j}) is not bonded in the molecule.")
        except Exception as e:
            print(str(e))   
            return {"dihedrals": True}

        
def get_bonds(connectivity, n_atoms):
    """Get a list of lists representing the bonds for each atom.

    Parameters
    ----------
    connectivity : list
        List of lists of length three containing the indices of the bonded atoms and the order of the bond.
    n_atoms : int
        Number of atoms in the molecule

    Returns
    -------
    list[list[int]]
        List of lists containing the indices of atoms connected to the atom represented by that index.
    """
    bonds = [[] for _ in range(n_atoms)]
    for i, j, _ in connectivity:
        if i == j:
            return []
        bonds[i].append(j)
        bonds[j].append(i)
    return bonds


def get_n_complexes(bonds, n_atoms):
    """Determine the number of components in the molecular object

    Parameters
    ----------
    bonds : list
        List of lists containing the indices of atoms connected to the atom represented by that index.
    n_atoms : int
        Number of atoms in the molecule(s)

    Returns
    -------
    int
        Number of separate graphs
    """
    graph = nx.Graph()
    for i in range(n_atoms):
        graph.add_node(i)

    for i, js in enumerate(bonds):
        for j in js:
            if not graph.has_edge(i, j):
                graph.add_edge(i, j)

    subgraphs = list(nx.connected_components(graph))
    return len(subgraphs)


def check_entry(entry, dataset_type):
    """Take a qcportal.external.scaffold entry and return keys representing tests that haven't passed.
    
    Note that LinearTorsionError are excluded from this test.

    Parameters
    ----------
    entry : dict
        QCFractal entry in dictionary form as exported from qcportal.external.scaffold
    dataset_type : str
        Type of dataset being assessed
        
    Returns
    -------
    errors : dict
        Keys indicate unpassed tests
    """
    errors = defaultdict()
    
    # Molecule Level
    key_mol = {
        "singlepoint": "molecule",
        "optimization": "initial_molecule",
        "torsiondrive": "initial_molecules",
    }.get(dataset_type)
    mol = entry[key_mol] if not isinstance(entry[key_mol], list) else entry[key_mol][0]
    if "molecular_charge" not in mol:
        errors["total_charge"] = True

    bonds = get_bonds(mol["connectivity"], int(len(mol["geometry"])/3)) if "connectivity" in mol else []
    if not bonds:
        errors["bonds"] = True
        warnings.warn("Bonds have not been defined, so the number of complexes and constraint definitions cannot be assessed.")
    else:
        try:
            n_complexes = get_n_complexes(bonds, int(len(mol["geometry"])/3))
            if n_complexes > 1:
                raise ValueError("More than one complex is present in the provided entry.")
        except Exception as e:
            print(str(e))
            errors["complex"] = True

    try:
        if "geometry" not in mol:
            raise ValueError("Molecular geometry is missing")
        else:
            coordinates = np.array(mol["geometry"]).reshape(-1, 3)
            distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=-1)
            distances[np.where(distances == 0)] = np.nan
            if np.min(distances) < 1.5: # 1.5 Bohr ~ 0.8 Angstroms
                raise ValueError("Bond distance less than 1.5 suggests units are not in Bohr.")
    except Exception as e:
        print(str(e))
        errors["coordinates"] = True
    
    # Entry Level
    ## cmiles
    if not [
        key for key in entry["attributes"].keys() 
        if glob.fnmatch.fnmatch(key, "canonical*hydrogen*smiles")
    ]:
        errors["cmiles"] = True

    ## Constraint formatting
    if "constraints" in entry["additional_keywords"]:
        errors.update(get_constraint_checks(entry["additional_keywords"]["constraints"], bonds=bonds))
        
    return errors