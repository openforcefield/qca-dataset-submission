from constructure.scaffolds import SCAFFOLDS
from utils import *
from collections import defaultdict

def process(input_file, data_name,store=True, filter=True, filter_type=2, check_n_rings=True, filter_ortho=True, isomeric=True):
    from collections import defaultdict
    from constructure.scaffolds import SCAFFOLDS
    from tqdm import tqdm

    frag_tot = defaultdict(list)

    for scaff_name in tqdm(SCAFFOLDS.keys()):
        frags = sdf_to_frag_smiles(sdf=input_file, scaffold=SCAFFOLDS[scaff_name], filter=filter, filter_type=filter_type, check_n_rings=check_n_rings, filter_ortho=filter_ortho, isomeric=isomeric)

        for rx, lst in frags.items():
            frag_tot[rx] = sorted(list(set(frag_tot[rx]+ frags[rx])))

    import itertools
    frag_list = [value for value in frag_tot.values()]
    frag_list = list(set(list(itertools.chain.from_iterable(frag_list))))

    if store:
        import pickle
        substituents = frag_list
        pname = data_name+ "_substituents.p"
        pickle.dump( substituents, open( pname , "wb" ) )

    groups = group_sub_by_type(frag_list)
    if store:
        pdfname = data_name+ "_substituents.pdf"
        draw_frag_pdf(groups, pdf_filename=pdfname)

    return frag_list, groups

# Low rank substituents: 
# 1. Hydrazines                                       : '*N1' w/ filter_hydrazine
# 2. Halogen attached to noncarbon (R1#101 and R1#102): '[!C;!c]-[F,Cl,Br]'
# 3. Pentavalent nitrogens (R1 #91)                   : '[N]=[N]#[N]'
# 4. Fluoren isotope (R1#65)                          : '[18F]'
# 5. Diff. halogen atoms attached to the same carbon  : 'C(F)[Cl,Br,I]', 'C(Cl)[Br,I]', 'C(Br)I'
def is_this_low_rank_substituent(smi):
    checklist = ['[!C;!c]-[F,Cl,Br,I]', '[N]=[N]#[N]', '[18F]','C(F)[Cl,Br,I]', 'C(Cl)[Br,I]', 'C(Br)I']
    match=False
    for target_smiles in checklist:
        if substructure_search(smi, target_smiles):
            match=True
    if search_hydrazine(smi):
        match=True

    return match

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

def filter_sub_list(sub_list, store=True):
    import pickle
    filtered = set()
    for canonical_smi in sub_list:
        if not is_this_low_rank_substituent(canonical_smi):
            if not substructure_search(canonical_smi):
                filtered.add(canonical_smi)

    frag_dict = group_sub_by_type(filtered)
    if store:
        draw_frag_pdf(frag_dict, scaffold=None, pdf_filename='substituents_filtered.pdf')
        pickle.dump(filtered, open('substituents_filtered.p','wb'))
    return filtered

def main():
    
    print('1. Roche set')
    roche_substituents, roche_groups = process(input_file='input-mol-sets/OpenFF_references.sdf', data_name='roche',store=True, isomeric=False)
    print('2. Coverage set')
    coverage_substituents, coverage_groups = process(input_file='input-mol-sets/chosen_supplemented.smi', data_name='coverage',store=True, isomeric=False)
    print('3. Pfizer discrepancy set')
    pfizer_substituents, pfizer_groups = process(input_file='input-mol-sets/PFE-OFF-100Frags.smi', data_name='pfizer',store=True, isomeric=False)
    print('4. Bayer set')
    bayer_substituents, bayer_groups = process(input_file='input-mol-sets/Bayer_no_Si.smi', data_name='bayer',store=True, isomeric=False)

    # combine into one list
    combined_list = list(set(roche_substituents+ coverage_substituents+ pfizer_substituents+bayer_substituents))
    print(f'Before removing duplicates: {len(roche_substituents+ coverage_substituents+ pfizer_substituents+bayer_substituents)}')
    print(f'Number of substituents: {len(combined_list)}')
    frag_dict = group_sub_by_type(combined_list)
    draw_frag_pdf(frag_dict, scaffold=None, pdf_filename='substituents_combined.pdf')
    
    import pickle
    pickle.dump(combined_list, open('substituents_combined.p','wb'))

    # remove low rank substituents ('substituents_filtered.p')
    filter_sub_list(combined_list)

if __name__=='__main__':
    main()