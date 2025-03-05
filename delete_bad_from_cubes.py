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



# Create a defaultdict with a default
# value of an empty list
field = 20
cubes_folder = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/slices/'%(field)
pruebas = '/home/data/alvaro/gns_gd/gns2/F%s/pruebas/'%(field)
clean = '/home/data/GNS/2021/H/%s/cleaned/'%(field)
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
    
    command = ['imcopy', clean + f'cube{k}.fits.gz[*,*,{ll}:{lh}]', pruebas + f'cube{k}_cl.fits']
    
    print(command)
    result = subprocess.run(command, check=True,cwd=pruebas )

    
    
    
    
    
    

