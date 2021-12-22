import copy
import json
import logging
from tempfile import NamedTemporaryFile

from openeye import oechem
from openff.qcsubmit.results.filters import (
    ResultRecordFilter,
)
from openff.toolkit.utils import UndefinedStereochemistryError


class UndefinedStereoFilter(ResultRecordFilter):
    def _filter_function(self, result, record, molecule) -> bool:

        has_stereochemistry = True

        molecule = copy.deepcopy(molecule)
        molecule._conformers = [molecule.conformers[0]]

        try:

            with NamedTemporaryFile(suffix=".sdf") as file:
                molecule.to_file(file.name, "SDF")
                molecule.from_file(file.name)

        except UndefinedStereochemistryError:
            has_stereochemistry = False

        return has_stereochemistry


def not_a_nitro_group(protomer):
    smarts = "[#8:1]=[#7:2]=[#8:3]"
    ss = oechem.OESubSearch(smarts)
    oechem.OEPrepareSearch(protomer, ss)
    if ss.SingleMatch(protomer):
        print("Yes, nitro group present")
        return False
    else:
        print("No, nitro group is not present")
        return True


def main():
    logging.basicConfig(level=logging.INFO)

    from openeye import oedepict

    multi = oedepict.OEMultiPageImageFile(oedepict.OEPageOrientation_Landscape,
                                          oedepict.OEPageSize_US_Letter)
    image = multi.NewPage()
    opts = oedepict.OE2DMolDisplayOptions()
    rows, cols = 3, 2
    grid = oedepict.OEImageGrid(image, rows, cols)
    grid.SetCellGap(20)
    grid.SetMargins(20)
    citer = grid.GetCells()

    with open("data-sets/1-2-0-opt-set.json", "r") as file:
        optimization_set = json.load(file)

    from openeye import oechem, oequacpac

    positive_charges = []
    net_charges = []
    negative_charges = []
    zero_charges = []
    proto_positive_charges = []
    proto_net_charges = []
    proto_negative_charges = []
    proto_zero_charges = []
    smiles = {}
    n = 0
    n_good = 0
    entries = list(optimization_set['entries'].values())[0]
    ofs = oechem.oemolostream('confs_with_different_protonation_states.smi')

    for entry in entries:
        # create oemol from canonical smiles
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, entry['cmiles'])
        oechem.OEAddExplicitHydrogens(mol, True, True)

        mol_natoms = mol.NumAtoms()
        # output the confs with protonation states

        smi = oechem.OEMolToSmiles(mol)
        mol_neg_charges = []
        mol_pos_charges = []

        if smi not in smiles:
            n = n + 1
            smiles[smi] = True
            mol_net_chg = 0
            for atom in mol.GetAtoms():
                atom_chg = atom.GetFormalCharge()
                mol_net_chg = mol_net_chg + atom_chg
                if atom_chg < 0:
                    mol_neg_charges.append(atom_chg)
                    negative_charges.append(atom_chg)
                elif atom_chg > 0:
                    mol_pos_charges.append(atom_chg)
                    positive_charges.append(atom_chg)
                elif atom_chg == 0:
                    zero_charges.append(atom_chg)
            net_charges.append(mol_net_chg)
            num_protomers = 0
            for protomer in oequacpac.OEGetReasonableProtomers(mol):
                num_protomers = num_protomers + 1

            atleast_one_reasonable_protomer_present = False
            for protomer in oequacpac.OEGetReasonableProtomers(mol):
                if not atleast_one_reasonable_protomer_present:
                    oechem.OEAddExplicitHydrogens(protomer)
                    prot_natoms = protomer.NumAtoms()
                    mol_proto_neg_charges = []
                    mol_proto_pos_charges = []

                    net_chg = 0
                    for atom in protomer.GetAtoms():
                        atom_chg = atom.GetFormalCharge()
                        net_chg = net_chg + atom_chg
                        if atom_chg < 0:
                            mol_proto_neg_charges.append(atom_chg)
                            proto_negative_charges.append(atom_chg)
                        elif atom_chg > 0:
                            mol_proto_pos_charges.append(atom_chg)
                            proto_positive_charges.append(atom_chg)
                        elif atom_chg == 0:
                            proto_zero_charges.append(atom_chg)

                    if num_protomers == 1:
                        if mol_natoms == prot_natoms:
                            print("good protonation state")
                            n_good = n_good + 1
                            proto_net_charges.append(net_chg)
                        else:
                            proto_net_charges.append(net_chg)
                            if not citer.IsValid():
                                # go to next page
                                image = multi.NewPage()
                                grid = oedepict.OEImageGrid(image, rows, cols)
                                grid.SetCellGap(20)
                                grid.SetMargins(20)
                                citer = grid.GetCells()

                            cell = citer.Target()
                            mol.SetTitle("Used in fit: natoms = " + str(mol_natoms) + ", net_chg = " + str(
                                mol_net_chg))
                            oedepict.OEPrepareDepiction(mol, True, False)
                            opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
                            disp = oedepict.OE2DMolDisplay(mol, opts)
                            oedepict.OERenderMolecule(cell, disp)
                            oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OERedPen))
                            citer.Next()

                            cell = citer.Target()
                            protomer.SetTitle(
                                "Reasonable protomer: natoms = " + str(prot_natoms) + ", net_chg = " + str(
                                    net_chg))
                            oedepict.OEPrepareDepiction(protomer, True, False)
                            opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
                            disp = oedepict.OE2DMolDisplay(protomer, opts)
                            oedepict.OERenderMolecule(cell, disp)
                            oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OERedPen))
                            citer.Next()
                            if not_a_nitro_group(protomer):
                                print("Writing protomer")
                                oechem.OEWriteMolecule(ofs, protomer)
                            atleast_one_reasonable_protomer_present = True
                    elif num_protomers > 1:
                        if mol_natoms == prot_natoms and not atleast_one_reasonable_protomer_present:
                            proto_net_charges.append(net_chg)
                            n_good = n_good + 1
                            print("atleast_one_reasonable_protomer_present")
                            atleast_one_reasonable_protomer_present = True

            if not atleast_one_reasonable_protomer_present:
                flag_added_one_protomer = False
                for protomer in oequacpac.OEGetReasonableProtomers(mol):
                    if not flag_added_one_protomer:
                        flag_added_one_protomer = True
                        oechem.OEAddExplicitHydrogens(protomer, True, True)

                        net_chg = 0
                        prot_natoms = protomer.NumAtoms()
                        for atom in protomer.GetAtoms():
                            atom_chg = atom.GetFormalCharge()
                            net_chg = net_chg + atom_chg
                        if mol_natoms != prot_natoms:
                            proto_net_charges.append(net_chg)
                            if not citer.IsValid():
                                # go to next page
                                image = multi.NewPage()
                                grid = oedepict.OEImageGrid(image, rows, cols)
                                grid.SetCellGap(20)
                                grid.SetMargins(20)
                                citer = grid.GetCells()

                            cell = citer.Target()
                            mol.SetTitle("Used in fit: natoms = " + str(mol_natoms) + ", net_chg = " + str(mol_net_chg))
                            oedepict.OEPrepareDepiction(mol, True, False)
                            opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
                            disp = oedepict.OE2DMolDisplay(mol, opts)
                            oedepict.OERenderMolecule(cell, disp)
                            oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OERedPen))
                            citer.Next()

                            cell = citer.Target()
                            protomer.SetTitle(
                                "Reasonable protomer: natoms = " + str(prot_natoms) + ", net_chg = " + str(
                                    net_chg))
                            oedepict.OEPrepareDepiction(protomer, True, False)
                            opts.SetDimensions(cell.GetWidth(), cell.GetHeight(), oedepict.OEScale_AutoScale)
                            disp = oedepict.OE2DMolDisplay(protomer, opts)
                            oedepict.OERenderMolecule(cell, disp)
                            oedepict.OEDrawBorder(cell, oedepict.OEPen(oedepict.OERedPen))
                            citer.Next()
                            if not_a_nitro_group(protomer):
                                print("Writing protomer")
                                oechem.OEWriteMolecule(ofs, protomer)

    import matplotlib.pyplot as plt

    plt.figure()
    plt.hist(negative_charges, alpha=0.5, rwidth=0.5, label='used in fit')
    plt.hist(proto_negative_charges, alpha=0.5, label='after enumerating protomers')
    plt.legend()
    plt.title("Gen 2 opt training set")
    plt.xlabel("Histogram of negative charges")
    plt.show()
    plt.savefig("hist_of_training_negative_charges.png", dpi=600)

    plt.figure()
    plt.hist(net_charges, alpha=0.5, label='used in fit', range=[-3, 3])
    plt.hist(proto_net_charges, alpha=0.5, label='after enumerating protomers', range=[-3, 3])
    plt.legend()
    plt.title("Gen 2 opt training set")
    plt.xlabel("Histogram of net charges")
    plt.show()
    plt.savefig("hist_of_training_net_charges.png", dpi=600)

    plt.figure()
    plt.hist(positive_charges, alpha=0.5, rwidth=0.5, label='used in fit')
    plt.hist(proto_positive_charges, alpha=0.5, label='after enumerating protomers')
    plt.legend()
    plt.title("Gen 2 opt training set")
    plt.xlabel("Histogram of positive charges")
    plt.show()
    plt.savefig("hist_of_training_positive_charges.png", dpi=600)

    oedepict.OEWriteMultiPageImage("mols_with_different_protonation_states.pdf", multi)

    # output: n:  1223 , n_good:  758
    print("n: ", n, ", n_good: ", n_good)


if __name__ == "__main__":
    main()
