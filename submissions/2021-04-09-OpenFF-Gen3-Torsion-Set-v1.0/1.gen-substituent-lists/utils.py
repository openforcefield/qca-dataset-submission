def sdf_to_frag_smiles(sdf, scaffold, filter=False, filter_type=2, check_n_rings=True, rm_duplicate=True, filter_ortho=True, isomeric=True):
    import os 
    from fragmenter import chemi
    from collections import defaultdict

    fragments_tot = defaultdict(list)

    # generate oemols from the input file (sdf or smi)
    oemols = chemi.file_to_oemols(sdf)
    for i, oemol in enumerate(oemols):
        fragments = oemol_to_frag_smiles(oemol, scaffold, filter=filter, filter_type=filter_type, check_n_rings=check_n_rings, filter_ortho=filter_ortho, isomeric=isomeric)
        for rx, smis in fragments.items():
            fragments_tot[rx] += list(set(smis)-set(fragments_tot[rx]))

    if rm_duplicate :
        import copy
        fragments_tot_copy = copy.deepcopy(fragments_tot)
        fragments_tot = defaultdict(list)
        fragments_tot[1] = sorted(fragments_tot_copy[1])
        for rx, lst in fragments_tot_copy.items():
            overlap = False
            overlapped_rx = [rx]
            for nrx, nlst in fragments_tot.items():
                if set(lst) == set(nlst): 
                    if nrx != rx:
                        overlap = True
                        overlapped_rx.append(nrx)
            if overlap:
                nu_name = ','.join([str(i) for i in set(overlapped_rx)])
                for rx_to_remove in overlapped_rx: 
                    if rx_to_remove in fragments_tot:
                        fragments_tot.pop(rx_to_remove)
                fragments_tot[nu_name] = sorted(lst)
            if not overlap: 
                fragments_tot[rx] = sorted(lst)

    return fragments_tot

def oemol_to_frag_smiles(oemol, scaffold, filter=True, filter_type=2, check_n_rings=True,filter_ortho=True, isomeric=True): 

    from collections import defaultdict
    from constructure.utilities.openeye import remove_duplicate_smiles

    # get bonds to break 
    bonds_to_break = get_atom_indices_from_oemol(oemol, scaffold)

    fragments = defaultdict(list)
    for idx, bonds_to_break in bonds_to_break.items(): 
        rx = idx + 1
        for bond_to_break in bonds_to_break:
            atom1, atom2 = bond_to_break
            atoms, bonds = get_atoms_bonds_oemol(oemol, atom1, atom2)
            atom_bond_set = to_atom_bond_set(oemol, atoms, bonds)
            fragment = atom_bond_set_to_mol(atom_bond_set, oemol) 
            # convert into smiles
            frag_smi = frag_to_smile(fragment, isomeric)  
            if frag_smi is None: 
                pass 
            else: 
                if filter: 
                    _complex, reason = check_frag_complexity(frag_smi, filter_type=filter_type, check_n_rings=check_n_rings,filter_ortho=filter_ortho)
                    if _complex: 
                        print(f'Skip complex substructure, {frag_smi}(reason: {reason})')
                    else:
                        fragments[rx].append(frag_smi)
                else:
                    fragments[rx].append(frag_smi)
        fragments[rx] = remove_duplicate_smiles(fragments[rx])
    return fragments

def get_atom_indices_from_oemol(oemol, scaffold): 
    from openeye import oechem
    import re
    from collections import defaultdict
    import cmiles
    # gen scaffold_oemol
    scaffold_smi = re.sub(r"\(\[R([1-9])+]\)", r"([*])", scaffold.smiles)
    scaff_oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(scaff_oemol, scaffold_smi)
    # check if the scaffold oemol has isomericity, if not generate a random stereoisomer for substructure search
    try: 
        explicit_h_smiles = cmiles.utils.mol_to_smiles(scaff_oemol, mapped=False)
    except: 
        print(f'No stereochemistry defined:{ cmiles.utils.mol_to_smiles(scaff_oemol, mapped=False, isomeric=False)}')
        print('-> Will enumerate a random stereoisomer for substructure search.')
        from fragmenter.states import _enumerate_stereoisomers
        scaff_oemol = _enumerate_stereoisomers(scaff_oemol, max_states=1)[0]
        explicit_h_smiles = cmiles.utils.mol_to_smiles(scaff_oemol, mapped=False)
    # regenerate scaffold_oemol with stereo information
    scaff_oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(scaff_oemol, explicit_h_smiles)

    # Get the bonds to break 
    bonds_to_break_scaff = {}
    count = 0
    for idx, atom in enumerate(scaff_oemol.GetAtoms()): 
        if atom.GetAtomicNum() == 0: # R group
            for bond in atom.GetBonds():
                nbor = bond.GetNbr(atom)
                bonds_to_break_scaff[count] = (atom.GetIdx(), nbor.GetIdx())
                count += 1

    # substructure search
    ss = oechem.OESubSearch(explicit_h_smiles)
    oechem.OEPrepareSearch(oemol, ss)

    # 
    bonds_to_break = defaultdict(list)
    for count, match in enumerate(ss.Match(oemol)):
        pattern_atoms = [ma.pattern.GetIdx() for ma in match.GetAtoms()]
        target_atoms  = [ma.target.GetIdx() for ma in match.GetAtoms()]
        # sort
        sorted_target_atoms = [x for _,x in sorted(zip(pattern_atoms, target_atoms))]

        # convert into atom indices
        for idx, bond in bonds_to_break_scaff.items(): 
            atom1, atom2 = bond
            atom1_idx = sorted_target_atoms[atom1]
            atom2_idx = sorted_target_atoms[atom2]
            #check if (atom1_idx, atom2_idx) is in-ring or not 
            oeatom1 = oemol.GetAtom(oechem.OEHasAtomIdx(atom1_idx))
            oeatom2 = oemol.GetAtom(oechem.OEHasAtomIdx(atom2_idx))
            bond = oemol.GetBond(oeatom1, oeatom2)
            if not bond:
                raise ValueError(f"{(atom1_idx, atom2_idx)} is a disconnected bond")

            if not bond.IsInRing():
                if (atom1_idx, atom2_idx) not in bonds_to_break[idx]: 
                    bonds_to_break[idx].append((atom1_idx, atom2_idx))
            # else: 
            #     print(f'# {(atom1_idx, atom2_idx)} is in-ring bond! ')
    return bonds_to_break

def get_atoms_bonds_oemol(oemol, atom1, atom2):
    atoms = [atom1]
    added = True
    while added: 
        added=False
        for atom in atoms:
            bonds_activated = [(bond.GetBgnIdx(), bond.GetEndIdx()) for bond in oemol.GetBonds() if atom in (bond.GetBgnIdx(), bond.GetEndIdx())]
            linked_atoms  = [ i if j==atom else j for i,j in bonds_activated]

            for latom in linked_atoms:
                if latom != atom2 and latom not in atoms:
                    atoms.append(latom)
                    added = True
    bonds = [(bond.GetBgnIdx(), bond.GetEndIdx()) for bond in oemol.GetBonds() if bond.GetBgnIdx() in atoms and bond.GetEndIdx() in atoms]                
    return atoms, bonds

def to_atom_bond_set(oemol, atoms, bonds):
    from openeye import oechem
    atom_bond_set = oechem.OEAtomBondSet()
    for a_idx in atoms:
        atom = oemol.GetAtom(oechem.OEHasAtomIdx(a_idx))
        atom_bond_set.AddAtom(atom)

    for b_tuple in bonds:
        a1 = oemol.GetAtom(oechem.OEHasAtomIdx(b_tuple[0]))
        a2 = oemol.GetAtom(oechem.OEHasAtomIdx(b_tuple[-1]))
        bond = oemol.GetBond(a1, a2)
        if not bond:
            raise ValueError("{} is a disconnected bond".format(b_tuple))
        atom_bond_set.AddBond(bond)

    return atom_bond_set

def atom_bond_set_to_mol(frag, oemol, adjust_hcount=False, RGroup=True):
    from openeye import oechem
    import warnings
    fragatompred = oechem.OEIsAtomMember(frag.GetAtoms())
    fragbondpred = oechem.OEIsBondMember(frag.GetBonds())

    fragment_oemol = oechem.OEMol()
    adjustHCount = adjust_hcount
    oechem.OESubsetMol(fragment_oemol, oemol, fragatompred, fragbondpred, adjustHCount, RGroup)

    oechem.OEAddExplicitHydrogens(fragment_oemol)
    # sanity check that all atoms are bonded
    for atom in fragment_oemol.GetAtoms():
        if not list(atom.GetBonds()):
            warnings.warn("Yikes!!! An atom that is not bonded to any other atom in the fragment. "
                            "You probably ran into a bug. Please report the input molecule to the issue tracker")
    # Perceive stereo and check that defined stereo did not change
    oechem.OEPerceiveChiral(fragment_oemol)
    oechem.OE3DToAtomStereo(fragment_oemol)
    oechem.OE3DToBondStereo(fragment_oemol) 

    return fragment_oemol

def frag_to_smile(fragment, isomeric=True): ##
    import cmiles
    from openeye import oechem

    mol_copy = oechem.OEMol(fragment)
    try: 
        explicit_h_smiles = cmiles.utils.mol_to_smiles(mol_copy, mapped=False)
    except: 
        print(f'No stereochemistry defined:{ cmiles.utils.mol_to_smiles(mol_copy, mapped=False, isomeric=False)}')
        print('-> Will enumerate a random stereoisomer for substructure search.')
        from fragmenter.states import _enumerate_stereoisomers
        mol_copy = _enumerate_stereoisomers(mol_copy, max_states=1)[0]
        from cmiles._cmiles_oe import has_stereo_defined
        if not has_stereo_defined(mol_copy): 
            print(f'{fragment} failed to gen stereo isomer!!')
            return None
        else: 
            explicit_h_smiles = cmiles.utils.mol_to_smiles(mol_copy, mapped=False)
    try: 
        cmiles_identifiers = cmiles.get_molecule_ids(explicit_h_smiles, toolkit='openeye') 
        if isomeric: 
            return cmiles_identifiers['canonical_isomeric_smiles']  
        else:
            return cmiles_identifiers['canonical_smiles']  
    except: 
        print(f'{fragment} failed to get inchi!!')
        return None

def check_frag_complexity(frag_smi, filter_type=2, check_n_rings=True, filter_ortho=True):
    from openeye import oechem
    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, frag_smi)    
    nrots = oechem.OECount(oemol, oechem.OEIsRotor())
    # print(f'{frag_smi}: {nrots}')
    num_components, component_membership = oechem.OEDetermineComponents(oemol)
    num_rings = oemol.NumBonds() - oemol.NumAtoms() + num_components
    if filter_type == 1: 
        if nrots >0: 
            return True, f' nrot: {nrots}'
        else:
            if check_n_rings and num_rings > 1: 
                    return True, f' nrings: {num_rings}'
            else:
                if filter_ortho and find_ortho_substituents(frag_smi):
                    return True, f'ortho substituent exists.'
                else: 
                    return False, 'pass'  
    elif filter_type == 2: 
        if nrots >1: 
            return True, f' nrot: {nrots}'
        else: # nrot = 0 or 1 
            if check_n_rings and num_rings > 1: 
                    return True, f' nrings: {num_rings}'
            # remain 1ring with nrot 0 or 1/ chain with nrot 0 or 1
            elif check_n_rings and num_rings == 1 and nrots ==1:
                return True, f' nrings: {num_rings}, nrots: {nrots}'
            else:
                if filter_ortho and find_ortho_substituents(frag_smi):
                    return True, f'ortho substituent exists.'
                else: 
                    return False, 'pass' 

def find_ortho_substituents(frag_smi):
    from openeye import oechem

    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, frag_smi)
    oechem.OEAddExplicitHydrogens(oemol)

    ortho = '[!#1]~!@[*;r]~;@[*;r]~!@[!#1]'
    ss = oechem.OESubSearch(ortho)
    oechem.OEPrepareSearch(oemol, ss)

    return ss.SingleMatch(oemol)

def draw_frag_pdf(frag_dict, scaffold=None, pdf_filename='fragments.pdf'): 
    from openeye import oechem, oedepict
    import re
    itf = oechem.OEInterface()
    PageByPage = True
    suppress_h = True
    rows = 7
    cols = 5
    ropts = oedepict.OEReportOptions(rows, cols)
    ropts.SetHeaderHeight(25)
    ropts.SetFooterHeight(25)
    ropts.SetCellGap(2)
    ropts.SetPageMargins(10)
    report = oedepict.OEReport(ropts)
    cellwidth, cellheight = report.GetCellWidth(), report.GetCellHeight()
    opts = oedepict.OE2DMolDisplayOptions(cellwidth, cellheight, oedepict.OEScale_Default * 0.5)
    opts.SetAromaticStyle(oedepict.OEAromaticStyle_Circle)
    pen = oedepict.OEPen(oechem.OEBlack, oechem.OEBlack, oedepict.OEFill_On, 1.0)
    opts.SetDefaultBondPen(pen)
    oedepict.OESetup2DMolDisplayOptions(opts, itf)
    if scaffold: 
        oemol = oechem.OEGraphMol()
        scaffold_smi = re.sub(r"\(\[R([1-9])+]\)", r"([*])", scaffold.smiles)
        oechem.OESmilesToMol(oemol, scaffold_smi)
        cell = report.NewCell()
        mol = oechem.OEMol(oemol)
        mol.SetTitle(f'{scaffold.smiles}')
        oedepict.OEPrepareDepiction(mol, False, suppress_h)
        disp = oedepict.OE2DMolDisplay(mol, opts)
        oedepict.OERenderMolecule(cell, disp)
        headerpen = oedepict.OEPen(oechem.OEWhite, oechem.OELightGrey, oedepict.OEFill_Off, 1.0)
        oedepict.OEDrawBorder(cell, headerpen)

    for rx, smis in frag_dict.items():
        for idx, smi in enumerate(smis):
            if smi != None: 

                # Create oemol
                oemol = oechem.OEGraphMol()
                oechem.OESmilesToMol(oemol, smi)

                # Render molecule
                cell = report.NewCell()
                mol = oechem.OEMol(oemol)
                mol.SetTitle(f'R{rx} #{idx+1}')
                oedepict.OEPrepareDepiction(mol, False, suppress_h)
                disp = oedepict.OE2DMolDisplay(mol, opts)

                oedepict.OERenderMolecule(cell, disp)

    oedepict.OEWriteReport(pdf_filename, report)


def group_sub_by_type(substructures):
    from collections import defaultdict, OrderedDict
    from openeye import oechem
    groups = defaultdict(list)
    for smi in substructures: 

        oemol = oechem.OEGraphMol()
        oechem.OESmilesToMol(oemol, smi) 

        num_components, component_membership = oechem.OEDetermineComponents(oemol)
        num_rings = oemol.NumBonds() - oemol.NumAtoms() + num_components
        if num_rings == 0:
            groups[1].append(smi) # aliphatic chains 
        elif num_rings > 0: 
            nraromsystems, parts = oechem.OEDetermineAromaticRingSystems(oemol)
            if nraromsystems == 0:
                groups[2].append(smi) # aliphatic rings
            elif nraromsystems > 0:
                for atom in oemol.GetAtoms():
                    if parts[atom.GetIdx()] == 1:
                        size = oechem.OEAtomGetSmallestRingSize(atom)
                if size == 6: 
                    groups[3].append(smi) # 6-membered aromatic rings
                elif size == 5: 
                    groups[4].append(smi) # 5-membered aromatic rings
                else: 
                    print(f'not included in any group, {smi}')
            else: 
                print(f'not included in any group, {smi}')


    for rx, lst in groups.items():
        groups[rx] = sorted(lst)
    groups = OrderedDict(sorted(groups.items()))
    print(f'1. Acyclic aliphatic: {len(groups[1])}, 2. aliphatic rings: {len(groups[2])}, 3. 6-membered aromatic rings:{len(groups[3])}, 4. 5-membered aromatic rings:{len(groups[4])}')
    return groups

def mol2_to_oemol(mol2_fnm):
    from openeye import oechem
    ifs = oechem.oemolistream(mol2_fnm)
    oemol = oechem.OEMol()
    oechem.OEReadMolecule(ifs, oemol)
    return oemol

def search_hydrazine(smi):
    if smi.startswith('*N1'):
        return True
    else:
        return False

def substructure_search(smi, target_smiles='[!C;!c]-[F,Cl,Br]'):
    from openeye import oechem
    ss = oechem.OESubSearch(target_smiles)

    oemol = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol, smi)
    oechem.OEAddExplicitHydrogens(oemol)
    oechem.OEPrepareSearch(oemol, ss)

    return ss.SingleMatch(oemol)

def find_molecule_to_blame(input_file,  filter=True, filter_type=2, check_n_rings=True, filter_ortho=True, isomeric=True):
    
    checklist = ['[!C;!c]-[F,Cl,Br]', '[N]=[N]#[N]', '[18F]','C(F)[Cl,Br,I]', 'C(Cl)[Br,I]', 'C(Br)I']
    from collections import defaultdict
    from constructure.scaffolds import SCAFFOLDS
    import os
    from fragmenter import chemi

    molecules= defaultdict(set)
    fragments = defaultdict(set)

    for scaff_name in SCAFFOLDS.keys():
        scaffold=SCAFFOLDS[scaff_name]
        oemols = chemi.file_to_oemols(input_file)
        for i, oemol in enumerate(oemols):
            fragments_ = oemol_to_frag_smiles(oemol, scaffold, filter=filter, filter_type=filter_type, check_n_rings=check_n_rings, filter_ortho=filter_ortho, isomeric=isomeric)
            fragments_list =  list(set(sum(fragments_.values(), [])))
            letsblame = False
            reason = None
            frag = None
            for fragment in fragments_list: 
                # check1
                for target_smiles in checklist: 
                    if substructure_search(fragment, target_smiles):
                        letsblame=True
                        reason = target_smiles
                        frag = fragment
                # check2
                if search_hydrazine(fragment):
                    letsblame=True
                    reason = '*N1'
                    frag = fragment
            if letsblame: 
                # print smiles of oemol? 
                molecules[reason].add(frag_to_smile(oemol, isomeric=False))
                fragments[reason].add(frag)
    return molecules, fragments 