"""
01_generateOptDS.py filters optimization data sets from QCA through fingerprint clustering.

The script can be ran by performing the following command:
    python 01_generateOptDS.py `QCA data set name`




The main functions in the script are:
load_DS : Loads the ds from .txt file or from QCA DS
paramUsage: Creates dictionary of parameters (keys) and molecules which use these parameters
selectDiverseMols: Returns a list of diverse molecules after clustering with DBSCAN based off graph similarity scores
clusterAdjustment (optional): Scans through various eps and min_samples to get desired number of clusters
makeJson: Creates .json file from final list of molecules for QCA input


By Jessica Maat (jmaat@uci.edu)

Contributors: David Mobley, Chris Bayly, Hyesu Jang, Lee-Ping Wang, Jeffrey Wagner, Josh Horton, Chaya Stern, John Chodera
"""

from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openeye.oegraphsim import *
from openeye.oechem import *
import numpy as np
import matplotlib.pyplot as plt
import qcportal as ptl
import fragmenter
import cmiles
import json
import logging
import tempfile
import gzip
import pickle
import random

def load_DS(mol_file):
    """
    Description -
    Loads in .txt file of smiles and returns a list of all of the smiles

    Input -
    mol_file: A .txt file of smiles for data set analysis

    Return -
    smiles: A list of smiles from mol_file

    """
    file = open(mol_file, 'r')
    text = file.readlines()
    smiles = [ line.split()[0] for line in text]
    return smiles



def load_DS_QCA(DSName):
    """
    Description-
    Loads in a QCA DS and return a list of smiles

    Input-
    DSName: Name of QCA optimization data set

    Return -
    smilesDS: List of smiles from the DS

    """
    client = ptl.FractalClient()
    ds = client.get_collection("OptimizationDataset", DSName)
    spec_name = ds.list_specifications().index[0]
    print(f"Loading TorsionDrive Scans from [ {DSName} ] spec [{spec_name}]")
    print(f"Found {len(ds.df)} data entries")
    # load torsiondrive record ids from the dataset
    map_record_id_entry_index = {}


    smilesDS=[]
    for entry_index in ds.df.index:
        data_entry = ds.get_entry(entry_index)
        smiles=data_entry.attributes['canonical_isomeric_explicit_hydrogen_smiles']
        #print(smiles)
        smilesDS.append(smiles)

    file1 = open("filterLog.txt","w")#write mode
    file1.write("This is the optimization data set:" + str(DSName)+  "\n")

    file1.write("This is the length of the original DS:" +str(len(smilesDS)) +"\n")
    file1.write("This is the length of the original DS without duplicates:" + str(len(set(smilesDS))) + "\n")
    file1.write("This is the list of the original smiles:" + str(smilesDS) + "\n")
    file1.close()


    return smilesDS


def loadSmilesSDF(fileName):
    """
    Takes in sdf file with smiles tag and compiles a list of smiles

    Input -
    fileName : the name of the sdf file e.g. "myfile.sdf"

    Returns -
    output : List of smiles strings
    """
    with open(fileName, 'r') as f:
        lines = f.readlines()
        output = []
        for index, line in enumerate(lines):
            if '<Canonical_Smiles>' in line:
                output.append(lines[index+1].strip('\n'))

    return output



def load_DS_QCA(DSName):
    """
    Description-
    Loads in a QCA DS and return a list of smiles

    Input-
    DSName: Name of QCA optimization data set

    Return -
    smilesDS: List of smiles from the DS

    """
    client = ptl.FractalClient()
    ds = client.get_collection("OptimizationDataset", DSName)
    spec_name = ds.list_specifications().index[0]
    print(f"Loading TorsionDrive Scans from [ {DSName} ] spec [{spec_name}]")
    print(f"Found {len(ds.df)} data entries")
    # load torsiondrive record ids from the dataset
    map_record_id_entry_index = {}


    smilesDS=[]
    for entry_index in ds.df.index:
        data_entry = ds.get_entry(entry_index)
        smiles=data_entry.attributes['canonical_isomeric_explicit_hydrogen_smiles']
        #print(smiles)
        smilesDS.append(smiles)

    return smilesDS



def paramUsage(smilesList, offxml):
    """
    Description -
    Reads in list of smiles and returns a dictionary of .offxml style parameters as keys
    and smiles of molecules as items

    Input -
    smilesList: A list of smiles
    offxml: The .offxml format force field that the parameters will be used with

    Return -
    anglebondDict: A dictionary of .offxml style parameters as keys and smiles of molecules that utilize
    parameters. The returned dictionary is only for bond and angle parameters e.g. 'a1', 'b2', etc.
    Note: The function can be modified to return a dictionary of torsion parameters.
    """

    # Initialize storage
    torsionDict = dict()
    anglebondDict = dict()


    # Let's label using our RC force field
    forcefield = ForceField(offxml)

    # Loop over smiles
    for smi in smilesList:

        # Create a simple molecule from SMILES and turn it into a topology.
        molecule = Molecule.from_smiles(smi, allow_undefined_stereo = True)
        topology = Topology.from_molecules([molecule])

        # Run the molecule labeling
        molecule_force_list = forcefield.label_molecules(topology)


        # Print out a formatted description of the parameters applied to this molecule
        for mol_idx, mol_forces in enumerate(molecule_force_list):
            for force_tag, force_dict in mol_forces.items():
                for (atom_indices, parameter) in force_dict.items():
                    pid = parameter.id

                    #create two seperate parameter usage dictionaries for (1) angle and bonds and (2) torsions
                    if "a" in pid or "b" in pid:
                        if not pid in anglebondDict:
                            anglebondDict[pid] = set()
                        anglebondDict[pid].add(smi)

                    #Uncomment this for torsion dictionary
                    #if "t" in pid:
                    #    if not pid in torsionDict:
                    #        torsionDict[pid] = set()
                    #    torsionDict[pid].add(smi)

    #Write out the angle and bond dictionary to "anglebond.p" file
    pickle.dump(anglebondDict, open( "anglebond.p", "wb" ) )

    return anglebondDict



def selectDiverseMols(paramDict, parameter, eps=0.5, min_samples=5):
    """
    Description -
    selectDiverseMols returns a limited set diverse of SMILES for a given parameter in a parameter dictionary.
    The methodology is first computed the graph similarity score and then using DBSCAN to create diverse clusters.

    Input -
    paramDict: A dictionary of .offxml style parameters as keys and smiles of molecules that utilize the
    parameters of the key
    parameter: An .offxml style parameter of interest e.g. 'a19'
    eps: Controls the maximum distance between two samples considered to be a neighbor of the other.
    Default is 0.5. This is based of previous experiments and scientific input from David Mobley & Chris Bayly.
    min_samples: Minimum number of samples near a compound for it to be considered a core point. Default = 5.

    Return -
    filtered_Smiles: A list of diverse and limited smiles for a given parameter

    """


    # Select fingerprint type
    # Chris Bayly/Krisztina Boda recommend tree fingerprints with default path length of 4 bonds
    fptype = OEFPType_MACCS166
    #fptype= OEFPType_Tree

    # Make fingerprints
    release_set_fp = []
    release_set_mols = []

    release_smi = paramDict[parameter]

    for smi in release_smi:
        mol = OEMol()
        OEParseSmiles(mol, smi)
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, fptype)
        release_set_fp.append(fp)
        release_set_mols.append(mol)



    import numpy as np
    N = len(release_smi)
    t_matrix = np.zeros( (N,N), np.float)

    comparison_fps = release_set_fp
    comparison_mols = release_set_mols

    for n in range(N):
        for m in range(N):
            if n != m: #Make diagonal elements be zero for more convenience
                t_matrix[n,m] = OETanimoto( comparison_fps[n], comparison_fps[m])

    #Make distance matrix run from 1 being far to 0 being close rather than the opposite (Tanimoto)
    dis_matrix = np.ones_like(t_matrix)-t_matrix
    #print(dis_matrix)
    from sklearn.cluster import DBSCAN
    from sklearn import metrics


    clustering = DBSCAN(eps, min_samples, metric="precomputed")

    # Fit clustering
    db = clustering.fit(dis_matrix)

    # Pull labels
    labels = db.labels_

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)



    #print('Estimated number of clusters: %d' % n_clusters_)
    #print('Total number of points: %d' % len(dis_matrix))
    #print('Estimated number of noise points: %d' % n_noise_)


    # Determine how many compounds are in each cluster.
    # Cluster "-1" is the "outliers"/noise points that are not in clusters.

    mols_by_cluster = {}
    cluster_nrs = set(labels)

    for label in cluster_nrs:
        mols_by_cluster[label] = []

        for (idx, thislabel) in enumerate(labels):
            if thislabel==label:
                mols_by_cluster[label].append(comparison_mols[idx])

        print("%d molecules in cluster %s" % (len(mols_by_cluster[label]), label))

    filtered_Smiles=[]


    for label, mols in mols_by_cluster.items():
        if label != -1:
            #sort molecules by size
            #sort_mols=sorted(mols, key=lambda x: x.NumAtoms())
            #select smallest molecule

            randomMolecule=random.choice(mols)
            molSmiles=OEMolToSmiles(randomMolecule)
            filtered_Smiles.append(molSmiles)

    if len(filtered_Smiles)==0:
        randomMol=random.choice(mols_by_cluster[-1])
        smilesRan=OEMolToSmiles(randomMol)
        filtered_Smiles.append(smilesRan)
        print(smilesRan)

    return filtered_Smiles



def clusterAdjustment(paramDict, parameter):
    """
    clusterAdjustment takes in parameter dictionary and parameter to scan through various eps and
    min_samples to determine the input values for DBSCAN to create cluster sizes of at least 5.

    Input -
    paramDict: Dictionary of parameters (keys) and molecules that use corresponding parameters
    parameter: Parameter of interest for clustering e.g. 'a1'

    Return -
    cluster: A list of selected filtered diverse molecules in smiles

    """
    if len(paramDict[parameter]) < 6:
        return paramDict[parameter]

    min_samples=5
    cluster=[]
    eps=1.0


    while len(cluster) <= 5:


        if eps < .11 and min_samples >=1:
            cluster=selectDiverseMols(paramDict, parameter, .5, min_samples)
            min_samples-=1
            if len(cluster) >= 3:
                return cluster

        if eps>0.11:
            cluster=selectDiverseMols(paramDict, parameter, eps, min_samples)
            eps-=.01


        if eps<.11 and min_samples==1:
            return paramDict[parameter]

    return cluster



def makeJson(smiles):
    """
    makeJson takes in a list of smiles strings and expands the tautomeric and isomeric state of the molecules
    and generates a .json file from these molecules.
    The functional also generates .smi files that record processed canonical smiles, duplicates, omega failures,
    cmiles failures, and skipped ions.

    input:
    smiles: List of smiles strings

    return:
    optSmiles: List of smiles that are used as optimization inputs in the .json file.
    """


    with tempfile.NamedTemporaryFile('w+', suffix='.smi') as tmp:
        #smiles = [smile+'\n' for smile in smiles]
        #tmp.writelines(smiles)
        for line in smiles:
            tmp.writelines(line+'\n')
        tmp.seek(0)
        temp_name = tmp.name
        print(tmp.name)
        oemols = fragmenter.chemi.file_to_oemols(temp_name)

        optimization_input = []
        processed_canonical_smiles = []
        skipped = []
        duplicates = [] # duplicate states
        omega_failures = []
        cmiles_failures = []

        # Write out SDF file of all conformations
        ofs = oechem.oemolostream('optimization_inputs.sdf')
    optimizationCount=0
    for mol in oemols:
        # Filter out single atom molecules
        if mol.GetMaxAtomIdx() == 1:
            skipped.append(cmiles.utils.mol_to_smiles(mol, mapped=False))
            continue

        # Expand protonation states and stereoisomers
        states = fragmenter.states.enumerate_states(mol, stereoisomers=False, tautomers=False)
        for s in states:
            # Some states have valences that rdkit does not accept.
            try:
                cmiles_ids = cmiles.get_molecule_ids(s)
            except:
                cmiles_failures.append(s)
                continue

            # Drop duplicates
            canonical_smiles = cmiles_ids['canonical_smiles']
            if canonical_smiles in processed_canonical_smiles:
                logging.info('Found duplicate canonical SMILES {}'.format(canonical_smiles))
                duplicates.append(canonical_smiles)
                continue
            else:
                processed_canonical_smiles.append(canonical_smiles)

            # Generate molecule using mapped SMILES
            mapped_smiles = cmiles_ids['canonical_isomeric_explicit_hydrogen_mapped_smiles']
            m = cmiles.utils.load_molecule(s)
            try:
                # Omega fails for some molecules.
                conformers = fragmenter.chemi.generate_conformers(m)
            except RuntimeError:
                logging.info('Omega failed to generate conformers for {}'.format(cmiles_ids['canonical_isomeric_smiles']))
                # Omega failed
                omega_failures.append(cmiles_ids['canonical_isomeric_smiles'])
                continue
            qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mapped_smiles) for conf in conformers.GetConfs()]
            optimization_input.append({'initial_molecules': qcschema_molecules,
                                       'cmiles_identifiers': cmiles_ids})
            optimizationCount+=len(qcschema_molecules)
            # Write to SDF
            oechem.OEWriteMolecule(ofs, conformers)

    with gzip.open('optimization_inputs.json.gz', 'w') as f:
        f.write(json.dumps(optimization_input, indent=2, sort_keys=True).encode('utf-8'))

    ofs.close()

    save_smiles(processed_canonical_smiles, 'optimization_inputs.smi')
    save_smiles(duplicates, 'duplicates.smi')
    save_smiles(omega_failures, 'omega_failures.smi')
    save_smiles(cmiles_failures, 'cmiles_failures.smi')
    save_smiles(skipped, 'skipped_ions.smi')
    print("Number of unique molecules optimized:" + str(len(oemols)))
    print("Final optimization count is:" + str(optimizationCount))

    file1 = open("finalCounts.txt","w")#write mode
    file1.write("Number of molecules optimized:" + str(len(oemols)) +'\n')
    file1.write("Final optimization count with expanded states is:" + str(optimizationCount) +'\n')
    file1.close()


    optSmiles=[]
    for mol in oemols:
        optSmiles.append(OEMolToSmiles(mol))

    return optSmiles


def save_smiles(smiles, filename):
    """Write smiles str to smi file"""
    with open(filename, 'w') as f:
        for smi in smiles:
            f.write(smi + '\n')


def make_param_histogram(param_id_counts, param_ids, letter, title):
    """
    Function that generates a histogram of parameter usage.

    input:
    param_id_counts: A list of number of times a parameter id was used
    param_ids: Corresponding list of the parameter ids e.g. ['a1', 'a2']
    letter: The corresponding .offxml letter code for the parameter type e.g. for angles, 'a'
    title: The title of the histogram plot.

    return:
    none
    """

    parm_ids = [ pid for pid in param_ids if pid[0]==letter]
    parm_ids.sort()
    counts_parms = [param_id_counts[parm_id] for parm_id in parm_ids]
    split = int(len(parm_ids)/2)

    indices = np.arange(len(parm_ids))
    fix, ax = plt.subplots(2,1,figsize=(16,5))
    ax[0].set_yscale('log', nonposy='clip')
    ax[1].set_yscale('log', nonposy='clip')

    rects2 = ax[0].bar(indices[0:split], counts_parms[0:split] )
    ax[0].set_ylabel('Count')
    ax[0].set_xticks( indices[0:split])
    ax[0].set_xticklabels( parm_ids[0:split], rotation=-60, ha='left')
    ax[0].set_xlim(indices[0], indices[split])
    plt.yscale('log',nonposy='clip')
    rects2 = ax[1].bar(indices[split:], counts_parms[split:])
    ax[1].set_ylabel('Count')
    ax[1].set_xticks( indices[split:])
    ax[1].set_xticklabels( parm_ids[split:], rotation=-60, ha='left')
    ax[1].set_xlim(indices[split], indices[-1]+1)

    ax[0].set_title(title)
    plt.savefig("finalParamUsage_hist.pdf")
    plt.show()


def getParamKeys(paramDict):
    """
    Gives a list of the parameter keys used in a dictionary of keys and molecules.
    The function also saves the used parameters in the data set to paramCoverage.txt

    input -
    paramDict: Dictionary with .offxml based parameters as keys and molecules as the items.

    return -
    params: List of parameters used
    """
    params=[]
    for key, item in paramDict.items():
        params.append(key)


    file1 = open("bondCoverage.txt","w")#write mode
    file2=open("angleCoverage.txt", "w")
    for par in params:
        if "b" in par:
           file1.write(str(par) + "\n")
        if "a" in par:
            file2.write(str(par) + "\n")
    file1.close()
    file2.close()

    return params


def main(DSName):
    smiles=loadSmilesSDF(DSName)
    dictofParams=paramUsage(smiles, 'openff_unconstrained-1.0.0-RC1.offxml')
    filtered_smiles=[]

    #filtered_set=set()
    #filename='anglebond.p'
    #infile = open(filename,'rb')
    #dictofParams = pickle.load(infile)
    #infile.close()
    filtLength=[]

    for key, mols in dictofParams.items():
        print(key)
        print(len(mols))
        finalmols=selectDiverseMols(dictofParams, key, eps=0.3, min_samples=3)
        print(len(finalmols))
        #filtered_set.add(set(finalmols))
        filtered_smiles.extend(finalmols)
        filtLength.append([len(key), len(mols), len(finalmols)])

    print(len(filtered_smiles))
    #print(len(filtered_set))
    #print(len(smiles))
    print(filtLength)
    print(filtered_smiles)
    print(type(filtered_smiles))
    smilesList=makeJson(filtered_smiles)

    getParamKeys(dictofParams)
    print(filtLength)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Read in string.')
    parser.add_argument('strings', type=str)
    args = parser.parse_args()
    print('This is the input data set:', args.strings)
    main(args.strings)

