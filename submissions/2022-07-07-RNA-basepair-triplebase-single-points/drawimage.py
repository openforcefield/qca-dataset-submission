import os, sys
import glob as glob
import numpy as np
import warnings
import pickle
from PIL import Image
from PyPDF2 import PdfFileMerger
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import rdBase
#print(rdBase.rdkitVersion)
warnings.simplefilter("ignore")



with open('rdmols_loopMOTIFS.pkl', 'rb') as db:
    rdmols_loopmotifs = pickle.load(db)    
with open('rdmols_basepairCATALOG.pkl', 'rb') as db:
    rdmols_basepair = pickle.load(db)
with open('rdmols_triplebaseDB.pkl', 'rb') as db:
    rdmols_triplebase = pickle.load(db)
rdmols = rdmols_loopmotifs + rdmols_basepair + rdmols_triplebase

print("# of entries")
print("loopMOTIFS: {}".format(len(rdmols_loopmotifs)))
print("basepairCATALOG: {}".format(len(rdmols_basepair)))
print("triplebaseDB: {}".format(len(rdmols_triplebase)))
print("-------------------")
print("total: {}".format(len(rdmols)))


# convert to 2D
for rdmol in rdmols:
    AllChem.Compute2DCoords(rdmol, clearConfs=True)


# Draw 2D molecule and save image
# RuntimeError is raised when all rdmols are save as image. Therefore, images will be saved every 1000 entries and merged them into a single file
# https://github.com/rdkit/rdkit/issues/2974
# https://gist.github.com/greglandrum/56a80e84676fc3a24250b821d0faac13
# https://stackoverflow.com/questions/65470233/attributeerror-image-object-has-no-attribute-save

x = 1000 
n = int(len(rdmols)/x) + 1
for i in range(n):
    startidx = i * x
    endidx = startidx + x
    filename = "tmp" + str(i) + ".pdf"
    
    print("save images for entries #{} to #{}".format(startidx, endidx))
    img = Draw.MolsToGridImage(rdmols[startidx:endidx], molsPerRow=8, returnPNG=False, maxMols=99999)
    img = img.convert('RGB')
    img.save(filename)

print("merge pdf")
pdf = glob.glob("tmp*.pdf")
pdf.sort()
merger = PdfFileMerger()
for p in pdf:
    merger.append(p)
    os.remove(p)
merger.write("dataset.pdf")
merger.close()
