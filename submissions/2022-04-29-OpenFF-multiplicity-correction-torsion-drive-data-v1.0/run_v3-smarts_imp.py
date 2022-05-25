from smarts_impropers_angle_update_ff_3 import *
import qcportal as ptl
import pickle

client = ptl.FractalClient()

smilesdict={'t130b' : ['C1=CC=C(C=C1)S(=O)(=O)NN', 'CS(=O)(=O)NN(CC)S(=O)(=O)C', 'CCOC1=CC=C(C=C1)S(=O)(=O)NN', 'CN(C)C1=CC=CC2=C1C=CC=C2S(=O)(=O)NN', 'CS(=O)(=O)NN(CC)S(=O)(=O)C'],
't132b' : ['CC(C)NC(=O)C1=CC=C(C=C1)CNNC', 'CCNNC', 'CCNNCS(=O)(=O)C', 'C1=CC=C(C=C1)CCNN', 'c1ccc(cc1)C(=O)CNN'],
't133':['C1C[N+](CN1)([N+]2(CCNC2)[O-])[O-]',
'C1C[NH+]([NH+]1[O-])[O-]',
'C1CC[N+](CC1)([N+]2(CCCCC2)[O-])[O-]', 'CC(=O)C(C)(C)[NH+]([N+](C)(CCNC(C)(C)C(=O)CN1C=CN=C1[N+](=O)[O-])[O-])[O-]', 'CC(=O)C(C)(C)NCC[N+](C)([NH+](C(C)(C)C(=O)C)[O-])[O-]'],
't133a':['CN(CCC(=O)[O-])[N+](C)(C)C', 'C[N+](C)(C)NCCC(=O)[O-]', 'C[N+](C)(C=C)NCCC(=O)[O-]', 'CN(C)[N+](C)(C)C', 'C[N+](C)(C)NCCS(=O)(=O)[O-]'],
't142b':['CNSN(C)C(=O)ON', 'CN1C(=O)C=CS1', 'CCCCCCCCN1C(=O)C=CS1', 'C1=CC=C2C(=C1)C(=O)NS2', 'C1CCC(CC1)NSC2=NC3=CC=CC=C3S2'],
't142d':['C1=N[S+](NC1=N)[O-]', 'CC(=O)O[S+](C)C', 'CCN=C1C(=N[S+](N1)[O-])N', 'C[S+](C)OS(=O)(=O)c1ccccc1', 'c1ccc(cc1)[S+](c2ccccc2)O'],
't142e':['CN1C(=O)c2cc(c(cc2[S+]1[O-])OC)OC', 'C1(=N[S+](NC1=N)[O-])N', 'CN=C1C(=NCCS)N[S+](N1)[O-]', 'CSCC(N=C1C(=N[S+](N1)[O-])N)O', 'c1cc(oc1)CN=C2C(=N[S+](N2)[O-])N'],
't142f':['c1ccc(cc1)[S+](c2ccccc2)S(=O)(=O)O', 'C[SH+]S(=O)(=O)c1ccccc1', 'COc1ccccc1S(=O)(=O)[S+](C)S(=O)(=O)C', 'CC[S+](CC)S(=O)(=O)c1ccc[nH]1', 'CC[S+](C)S(=O)(=O)O'],
't142c':['CC(C(=O)OC)SP(=S)(OC)OC', 'COP(=S)(OC)SCN1C(=O)c2ccccc2C1=O', 'CCOP(=S)(OCC)S', 'CCSP(=S)(OC(C)C)OC(C)C', 'COP(=S)(OC)Sc1ccccc1'],
't143d':['c1ccc(cc1)C2=NS(=O)ON2', 'C1=NS(=O)ON1', 'Cc1ccs(=O)n1', 'CC(C(=O)C1=NS(=O)ON1)O', 'CCc1cs(=O)nc1N'],
't143e':['CN(C)S(=O)N(C)c1ccccc1', 'CCOS(=O)NC', 'CCOS(=O)Nc1ccccc1', 'CN(C)S(=O)N(C)C', 'NS(=O)O'],
't143f':['C[N+]1(CCCC1)S(=O)O', 'C[N+](C)(C)S(=O)O', 'C[N+](C)(c1ccccc1)S(=O)O', '[NH3+]S(=O)O', 'CC(=O)[NH2+]S(=O)O'],
't122d':['CC1=CC=C(C=C1)S(=O)C#C', 'CC(C)(C)OC#CS(=O)C1=CC=CC=C1', 'C1=CC=C(C=C1)C2=CN=[C-]S2=O', 'CCCC#CS(=O)C1=CC=C(C=C1)C', 'c1ccc(cc1)S(=O)C#N'],
't122e':['c1ccc(cc1)S(=O)c2ccccc2', 'c1ccc(cc1)S(=O)c2ccc3c(c2)[nH]c(n3)N', 'c1ccc2c(c1)Sc3ccccc3S2=O', 'Cc1ccc2c(c1)N(C(=O)c3ccccc3S2=O)C', 'c1cc(cc(c1)S(=O)c2cc(cc(c2)N)N)N'],
't122a': ['CC#CS(=O)(=O)C1=CC=CC=C1', 'C1=CC=C(C=C1)S(=O)(=O)C#N', 'CC1=CC=C(C=C1)S(=O)(=O)C#CC2=CC=CC=C2', 'CC#CS(=O)(=O)c1cc(cc(c1)N)N', 'c1c(cc(cc1O)S(=O)(=O)C#N)N'],
't116b':['C=C[S+]1CCCC1', 'C[S+](C)C=CC(=O)[O-]', 'CCC(C=C[S+](C)C)N', 'C[S+](C)C=CC(=O)Nc1ccccc1', 'C=C[S+]1CC(C(C1)N)N'],
't116c':['CN1C2C[S+]3CCCC3C2N(C1=O)C', 'C[S+](C)CCC(C(=O)O)N', 'C[S+](C)CCC(=O)[O-]', 'C[S+](C)CCC(=O)Nc1ccccc1', 'C[S+](C)CCc1ccccc1'],
't133':['C1C[N+](CN1)([N+]2(CCNC2)[O-])[O-]', 'C1C[NH+]([NH+]1[O-])[O-]', 'C1CC[N+](CC1)([N+]2(CCCCC2)[O-])[O-]', 'CC[N+](C)([NH+](C(C)(C)C(=O)C)[O-])[O-]', 'CC[NH+]([N+](C)(CC)[O-])[O-]'],
't33a':['CN(CCC(=O)[O-])[N+](C)(C)C', 'C[N+](C)(C)NCCC(=O)[O-]', 'C[N+](C)(C=C)NCCC(=O)[O-]', 'CN(C)[N+](C)(C)C', 'C[N+](C)(C)NCCS(=O)(=O)[O-]'],
't74a':['CC(=O)[N-]NC', 'CCC(=O)[N-]C(=O)C', 'CCC(=O)[N-]C', 'CS(=O)(=O)[N-]C(=O)c1ccc(cc1)N', 'CC[N-]C(=O)c1ccncc1'],
't157':['NS(=O)(=O)O', 'C(CS(=O)(=O)O)N', 'CS(=O)(=O)OCCCCOS(=O)(=O)C', 'OS(=O)(=O)O', 'CC1=CC(=O)NS(=O)(=O)O1'],
't157a': ['OS(=O)O', 'C(CS(=O)O)N', 'C(=N)(N)S(=O)O', 'CCOS(=O)OCC', 'C1COS(=O)O1']
}


smiles=[]

for key, item in smilesdict.items():
    smiles.extend(item)


get_gopt_matching_improper(smiles, 'openff-2.0.0-multiplicity-corrections-v8.offxml')





