#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:18:34 2024

@author: amartinez
"""

import astropy.io.fits as fits
import numpy as np
import os
import sys
from astropy.table import Table
from astropy.wcs.utils import fit_wcs_from_points
from astropy.wcs import WCS
from astropy.io import fits
import astroalign as aa
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs.utils import fit_wcs_from_points
import shutil
import gzip
import subprocess
import glob
# %
# %%
field = 20


folder = '/home/data/raw/GNS_2/H/Field/%s/'%(field)
# cubes_aligned = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/'%(field)
cubes_folder = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/slices/'%(field)

pruebas = '/home/data/alvaro/gns_gd/gns2/F%s/pruebas/'%(field)
sf_folder = '/home/data/GNS/2021/H/%s/data/'%(field)
clean = '/home/data/GNS/2021/H/%s/cleaned/'%(field)
VVV_fol = '/home/data/VVV/'
vvv_f = '/home/data/alvaro/gns_gd/gns2/VVV_fields/'
ims = '/home/data/GNS/2021/H/20/ims/'
tmp = '/home/data/alvaro/gns_gd/gns2/F20/tmp/'

# clean = pruebas
# %%
###########IMPORTANT################
# list of raw images and raw images are in /home/data/raw/GNS_2/H/Field/20
# list of cleand cubes is in '/home/data/GNS/2021/H/20/cleaned'
# list of starfinder stars is in '/home/data/GNS/2021/H/20/data'
# mask for all the chip is un '/home/data/GNS/2021/H/20/ims/mask.fits'
####################################

lista = open(folder + 'list.txt', 'r')
lines = len(lista.readlines())# This is closing lista. I have to reopen it
lista = open(folder + 'list.txt', 'r')
count = 0

with open(cubes_folder + '%s_cubes_slices_ext.txt'%(field),'w') as fil:
    fil.write('Cube_id #slices sl ext\n')

# lista = [1,2]
for i, li in  enumerate(lista):

    # print(li)
    name =  f'cube{i+1}.fits'
    command = ['gzip','-dk',clean +  f'{name}.gz']
    result = subprocess.run(command, check=True, cwd = tmp)

    os.replace(clean + name, tmp + name)
    
    
    
    command = ['fitsheader', name, '-k', 'NAXIS3']
    result = subprocess.run(command, check=True, capture_output=True, text=True, cwd = tmp)
    # Extract the value from the output
    output_line = result.stdout.strip()
    value = output_line.split('=')[1].split('/')[0].strip()
    
    # Convert the value to an integer (if needed)
    slices = int(value)
    
   
    print(30*'*' +f'\nSlices ={slices}\n' + 30*'*')
    command = ['missfits', name, '-OUTFILE_TYPE', 'SLICE', '-SAVE_TYPE','replace','-SLICE_SUFFIX','.%04d.fits', '-VERBOSE_TYPE','FULL'   ]
    result = subprocess.run(command, check=True,cwd=tmp)
    
    
    for sl in range(1,slices+1):
        count += 1 
        with open(cubes_folder + '%s_cubes_slices_ext.txt'%(field),'a') as fil:
            fil.write('%s %s %s %s\n'%(i+1,slices,sl,count))
        for chip in range(1,5):    
            
            name_c =  f'cubes_f{field}_c{chip}.{count:04d}.fits' 
            if chip ==1:
                command = ['fitscopy', f'cube{i+1}.{sl:04d}.fits[1:2048,1:2048]',name_c]
                result = subprocess.run(command, check=True,cwd=tmp)
                
            if chip ==2:
                command = ['fitscopy', f'cube{i+1}.{sl:04d}.fits[2049:4096,1:2048]', name_c]
                
                result = subprocess.run(command, check=True,cwd=tmp)
            if chip ==3:
                command = ['fitscopy', f'cube{i+1}.{sl:04d}.fits[2049:4096,2049:4096]', name_c]
                result = subprocess.run(command, check=True,cwd=tmp)
            if chip ==4:
                command = ['fitscopy', f'cube{i+1}.{sl:04d}.fits[1:2048,2049:4096]', name_c]
                result = subprocess.run(command, check=True,cwd=tmp)
        
        f_rm = os.path.join(tmp, f'cube{i+1}.{sl:04d}.fits')            
        os.remove(f_rm)
        print(f'Removed {f_rm}')
    
   
    
    
# 
for chip in range(1,5):
  
    # command = ['missfits', f'cubes_f{field}_c{chip}', '-outfile_type', 'cube', '-SAVE_TYPE', 'REPLACE', '-SLICE_SUFFIX','.%04d.fits']  
    # command = ['missfits', f'cubes_f{field}_c{chip}', '-c', 'default_2.missfits']  

    command = ['missfits', f'cubes_f{field}_c{chip}', '-outfile_type', 'cube','-SLICE_SUFFIX','.%04d.fits' ,'-SAVE_TYPE', 'REPLACE']  
    result = subprocess.run(command, check=True,cwd=tmp)
    
    pattern = os.path.join(tmp, f'cubes_f{field}_c{chip}.0*')
    
    # Find and remove matching files
    for file_path in glob.glob(pattern):
        try:
            os.remove(file_path)
            print(f"Removed: {file_path}")
        except OSError as e:
            print(f"Error removing {file_path}: {e}")





        





        
