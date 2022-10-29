import os, sys
import glob as glob
import numpy as np
import warnings
import pickle
from PIL import Image
from PyPDF2 import PdfFileMerger
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole, DrawingOptions
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import rdBase
#print(rdBase.rdkitVersion)
warnings.simplefilter("ignore")



with open('rdmols_loopMOTIFS.pkl', 'rb') as db:
    rdmols = pickle.load(db)    

print("# of entries")
print("total: {}".format(len(rdmols)))


# convert to 2D
mylist = []
_rdmols = []
for rdmol in rdmols:
    inchi = Chem.MolToInchi(rdmol)
    if inchi not in mylist:
        _rdmols.append(rdmol)
        AllChem.Compute2DCoords(_rdmols[-1], clearConfs=True)
        mylist.append(inchi)


# Draw 2D molecule and save image
# RuntimeError is raised when all rdmols are save as image. Therefore, images will be saved every 1000 entries and merged them into a single file
# https://github.com/rdkit/rdkit/issues/2974
# https://gist.github.com/greglandrum/56a80e84676fc3a24250b821d0faac13
# https://stackoverflow.com/questions/65470233/attributeerror-image-object-has-no-attribute-save

IPythonConsole.drawOptions.minFontSize=12
img = Draw.MolsToGridImage(_rdmols, molsPerRow=5, subImgSize=(500,500), returnPNG=False, maxMols=99999)
img = img.convert('RGB')
img.save("dataset.pdf")
