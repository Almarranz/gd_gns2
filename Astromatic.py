#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 12:54:33 2025

@author: amartinez
"""


import os
import sys
import shutil
import subprocess
import time
import glob

# %%

field = 20


folder = '/home/data/raw/GNS_2/H/Field/%s/'%(field)

pruebas = '/home/data/alvaro/gns_gd/gns2/F%s/pruebas/'%(field)
sf_folder = '/home/data/GNS/2021/H/%s/data/'%(field)
clean = '/home/data/GNS/2021/H/%s/cleaned/'%(field)
VVV_fol = '/home/data/VVV/'
ims = '/home/data/GNS/2021/H/20/ims/'


           
sex_folder = '/home/data/alvaro/gns_gd/gns2/F%s/sextractor/'%(field)
scamp_folder = '/home/data/alvaro/gns_gd/gns2/F%s/scamp/'%(field)
SWarp_folder = '/home/data/alvaro/gns_gd/gns2/F%s/SWarp/'%(field)
cubes_aligned = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/'%(field)
scripts = '/home/data/alvaro/gns_gd/gns2/scripts/'
out_folder = SWarp_folder + 'outputs/'

# %%
# %%
ch_range = [1,5]
#SOURCE-EXTRACTOR
t0_sex = time.time()
for chip in range(ch_range[0],ch_range[1]):
    command = ['source-extractor', cubes_aligned + '%s_image_c%s.fits'%(field,chip), 
               '-c', scripts + 'default_c.sex', '-CATALOG_NAME',f'{sex_folder}chip{chip}/%s_image_c%s.cat'%(field, chip),
               '-CHECKIMAGE_NAME',f'{sex_folder}chip{chip}/objects_%s_c%s.fits'%(field, chip)]
    
    try:
        # Run the command
        
        # result = subprocess.run(command, cwd=f'{sex_folder}chip{chip}/',check=True, text=True, capture_output=True)
        # result = subprocess.run(command, cwd=f'{sex_folder}chip{chip}/',check=True)
        result = subprocess.run(command, cwd= scripts,check=True)
        
        # Print standard output and error
        print("Command Output:")
        print(result.stdout)
        print("Command Error (if any):")
        print(result.stderr)
    
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Error: {e}")
        print(f"Standard Output: {e.stdout}")
        print(f"Standard Error: {e.stderr}")
t1_sex = time.time()

t_sex = t1_sex-t0_sex
print(30*'_' + '\nDone with SExtractor\nIt took %.0f hours\n'%(t_sex/3600) + 30*'_')

# %%
#SCAMP
# =============================================================================
# ASTREF_CATALOG         2MASS           # NONE,FILE,USNO-A2,USNO-B1,GSC-2.3,
#                                        # TYCHO-2,UCAC-4,URAT-1,NOMAD-1,PPMX,
#                                        # CMC-15,2MASS,DENIS-3,SDSS-R9,SDSS-R12,
#                                        # IGSL,GAIA-DR1,GAIA-DR2,GAIA-EDR3,
#                                        # PANSTARRS-1, or ALLWISE
# =============================================================================
t0_sca = time.time()
# for chip in range(ch_range[0],ch_range[1]):
for chip in range(2,5):
    command = ['scamp', sex_folder + 'chip%s/%s_image_c%s.cat'%(chip, field,chip), 
                '-c', 'scamp_c.conf', '-HEADER_NAME', f'{scamp_folder}chip{chip}/%s_image_c%s.head'%(field, chip)
                ,'-FULLOUTCAT_NAME',f'{scamp_folder}chip{chip}/%s_full_c%s.cat'%(field, chip),
                '-WRITE_XML','Y', '-XML_NAME', f'{scamp_folder}chip{chip}/scamp_f%s_c%s.xml'%(field, chip)]
    destination_folder = f'{scamp_folder}chip{chip}/'
    
    # command = ['scamp', sex_folder + 'chip%s/%s_image_c%s.cat'%(chip, field,chip), 
    #             '-c', 'scamp_c.conf', '-HEADER_NAME', f'{scamp_folder}chip{chip}__EGDR3/%s_image_c%s.head'%(field, chip)
    #             ,'-FULLOUTCAT_NAME',f'{scamp_folder}chip{chip}/%s_full_c%s.cat'%(field, chip),
    #             '-ASTREF_CATALOG','GAIA-EDR3']
    # destination_folder = f'{scamp_folder}chip{chip}_EGDR3/'
    
    try:
        # Run the command
        
        result = subprocess.run(command, cwd=scripts,check=True)
        # Print standard output and error
        print("Command Output:")
        print(result.stdout)
        print("Command Error (if any):")
        print(result.stderr)
    
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Error: {e}")
        print(f"Standard Output: {e.stdout}")
        print(f"Standard Error: {e.stderr}")
    fits_files = glob.glob(os.path.join(scripts, '*.pdf'))  # Find all .pdf files
    
    
    for fits_file in fits_files:
        try:
            # Construct destination file path
            if not os.path.exists(destination_folder):
                os.makedirs(destination_folder)
            dest_file = os.path.join(destination_folder, os.path.basename(fits_file))
            
            # Remove the file if it already exists in the destination folder
            if os.path.exists(dest_file):
                os.remove(dest_file)
                print(f"Existing file removed: {dest_file}")
            
            # Move the file
            shutil.move(fits_file, destination_folder)
            print(f"Moved: {fits_file} -> {dest_file}")
        except Exception as e:
            print(f"Error moving {fits_file}: {e}")
t1_sca = time.time()
t_sca = t1_sca-t0_sca
print(30*'_' + '\nDone with Scamp\nIt took %.0f min\n'%(t_sca/60) + 30*'_')

# %%
#SWARP
t0_swa = time.time()
# for chip in range(ch_range[0],ch_range[1]):
for chip in range(1,2):
    
    # command = ['SWarp', cubes_aligned+ '%s_image_c%s.fits' %(field, chip), 
    #            '-c', 'default_c.swarp', '-HEADER_NAME',scamp_folder + 'chip%s/%s_image_c%s.head'%(chip,field, chip),
    #            '-WEIGHT_IMAGE',cubes_aligned + '%s_mask_c%s.fits'%(field, chip),
    #            '-RESAMPLE_DIR', out_folder + 'chip%s/'%(chip)]
    command = ['SWarp', cubes_aligned+ '%s_image_c%s.fits' %(field, chip), 
               '-c', 'default_c.swarp', '-HEADER_NAME',scamp_folder + 'chip%s/%s_image_c%s.head'%(chip,field, chip),
               '-WEIGHT_IMAGE',cubes_aligned + '%s_mask_c%s.fits'%(field, chip),
               '-RESAMPLE_DIR', out_folder + 'chip%s/'%(chip),
               '-RESAMPLE','N',
               '-COMBINE','Y']
    
    try:
        # Run the command
        
        # result = subprocess.run(command, cwd=f'{SWarp_folder}chip{chip}/',check=True)
        result = subprocess.run(command, cwd= scripts ,check=True)
        # Print standard output and error
        print("Command Output:")
        print(result.stdout)
        print("Command Error (if any):")
        print(result.stderr)
    
    except subprocess.CalledProcessError as e:
        # Handle errors
        print(f"Error: {e}")
        print(f"Standard Output: {e.stdout}")
        print(f"Standard Error: {e.stderr}")
t1_swa = time.time()
t_swa = t1_swa-t0_swa
print(30*'_' + '\nDone with SWarp\nIt took %.0f hours\n'%(t_swa/3600) + 30*'_')