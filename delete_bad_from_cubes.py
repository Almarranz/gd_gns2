# %%
#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:38:03 2024

@author: alvaro
"""
from astropy.io import fits
import os
from collections import defaultdict
import numpy as np
from astropy.table import Table
import subprocess
import sys


# Create a defaultdict with a default
# value of an empty list
field = 20
cubes_folder = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/slices/'%(field)
pruebas = '/home/data/alvaro/gns_gd/gns2/F%s/pruebas/'%(field)
# clean = '/home/data/GNS/2021/H/%s/cleaned/'%(field)
clean_tmp = '/home/data/alvaro/gns_gd/gns2/F20/tmp/cleaned_tmp/'
sl = Table.read(cubes_folder + '20_cubes_slices_ext.txt', format = 'ascii')

print(sl.columns)

dic_bad =defaultdict(list)
dic_sl = defaultdict(list)
bad = [42,	47,	82,	222,	223,	241,	277,	386,	387,]


for b in bad:
    ind = np.where(b == sl['ext'])[0][0]
    # print(ind)
    dic_bad[sl[ind]['Cube_id']].append(sl[ind]['sl'])
    dic_sl[sl[ind]['Cube_id']].append((sl[ind]['#slices']))             





for k in list(dic_bad.keys()):
    # a_sl = sl['Cube_id'] == k
    a_sl = sl['#slices'][sl['Cube_id'] == k][0]
    sls = np.arange(1,a_sl+1)
    print(dic_bad[k])
    print(sls)
   
    mask = np.logical_not(np.isin(sls,dic_bad[k]))
    sls_g = sls[mask]
    print(sls_g)
    ll = min(sls_g)
    lh = max(sls_g)
    print(ll,lh)
    
    command = ['imcopy', clean_tmp + f'cube{k}.fits.gz[*,*,{ll}:{lh}]', clean_tmp + f'cube{k}.fits']
    print(command)
    result = subprocess.run(command, check=True,cwd=pruebas )
    
    
    cl_rm = os.path.join(clean_tmp,  f'cube{k}.fits.gz')
    os.remove(cl_rm)
    
    
command = 'gzip -d *.gz'
result = subprocess.run(command, check=True,cwd=clean_tmp, shell = True )

# Check the result (optional)
if result.returncode == 0:
    print("Successfully unzipped all .gz files in the directory.")
else:
    print("Failed to unzip .gz files.")    
    
    

