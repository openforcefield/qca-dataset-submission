#!/usr/bin/env python

"""
Retrieve current snapshot of SureChEMBL

"""

# Read contents of SureChEMBL directory
from urllib.request import urlopen
base_url = 'ftp://ftp.ebi.ac.uk/pub/databases/chembl/SureChEMBL/data/'
contents = urlopen(base_url).read().decode('utf-8')

# Identify all compressed SMILES files to retrieve
import re
pattern = re.compile('SureChEMBL_\d+_\d+.txt.gz') #the pattern actually creates duplicates in the list
filelist = pattern.findall(contents)
print(filelist)

# Retrieve files
import os
destination_directory = 'SureChEMBL'
if not os.path.exists(destination_directory):
    os.makedirs(destination_directory)
for filename in filelist:
    url = base_url + filename
    print(f'Retrieving {url}...')
    remotefile = urlopen(url)
    local_filename = os.path.join(destination_directory, filename)
    localfile = open(local_filename,'wb')
    localfile.write(remotefile.read())
    localfile.close()
    remotefile.close()
