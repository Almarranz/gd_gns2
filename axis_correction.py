# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Thu Oct 24 11:35:06 2024

@author: amartinez
"""

# %%
# Modify the size of the axes from the output of SWarp to make them equal
import numpy as np
from astropy.io import fits
import os
from astropy.table import Table
import sys
# pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'

field = 20
target_size = 2048
# %%
for chip in range(3, 5):
    folder = '/home/data/alvaro/gns_test/F%s/SWarp/outputs/chip%s/'%(field, chip)
    # folder = pruebas
    fits_files = [f for f in sorted(os.listdir(folder)) if f.endswith('.fits') and f.startswith('%s_image_c%s' % (field, chip)) or f.endswith('weight.fits') and f.startswith('%s_image_c%s' % (field, chip))]
    
    for nf, f_file in enumerate(fits_files):
        hdul = fits.open(folder + f_file)
        image_data = hdul[0].data
        ejes = [hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2']]
        print(f'Working with {f_file}')
        print(f'Original axes: {ejes}')
        
        arget_size = 2048
        axis1, axis2 = ejes

        # Calculate crop/padding amounts for each axis
        crop_pad_axis1 = (axis1 - target_size) // 2
        crop_pad_axis2 = (axis2 - target_size) // 2

        # Process for axis 1 if not equal to target size
        if axis1 != target_size:
            if axis1 > target_size:  # Crop
                start_x = crop_pad_axis1
                end_x = start_x + target_size
                image_data = image_data[:, start_x:end_x]
                hdul[0].header['CRPIX1'] -= start_x
            else:  # Pad
                pad_x = abs(crop_pad_axis1)
                image_data = np.pad(image_data, ((0, 0), (pad_x, pad_x)), mode='constant', constant_values=0)
                hdul[0].header['CRPIX1'] += pad_x

        # Process for axis 2 if not equal to target size
        if axis2 != target_size:
            if axis2 > target_size:  # Crop
                start_y = crop_pad_axis2
                end_y = start_y + target_size
                image_data = image_data[start_y:end_y, :]
                hdul[0].header['CRPIX2'] -= start_y
            else:  # Pad
                pad_y = abs(crop_pad_axis2)
                image_data = np.pad(image_data, ((pad_y, pad_y), (0, 0)), mode='constant', constant_values=0)
                hdul[0].header['CRPIX2'] += pad_y

        # Update header axes and save to file
        hdul[0].header['NAXIS1'] = target_size
        hdul[0].header['NAXIS2'] = target_size
        # hdul[0].data = image_data
        if f_file.endswith('weight.fits'):
            image_data = np.where(image_data ==0,0,1)
            hdul[0].data = image_data
        else:
            hdul[0].data = image_data
        hdul.writeto(folder + f'PADDED_{f_file}', overwrite=True)
        
        print(f'New axes: {hdul[0].header["NAXIS1"]}, {hdul[0].header["NAXIS2"]}')
        print(30 * '_')

# %%       
cubes_list = '/home/data/alvaro/gns_test/F%s/cubes_aligned/'%(field)
sl = Table.read(cubes_list + '%s_cubes_and_slices.txt'%(field), format = 'ascii')

ax_sz = target_size 

# chip = 1
for chip in range(3,5):
    folder = '/home/data/alvaro/gns_test/F%s/SWarp/outputs/chip%s/'%(field, chip)
    gd_folder = '/home/data/GNS/2021/H/%s/cubes_gd/chip%s/'%(field, chip)
    
    # for n,fi in enumerate(fits_files):
    idex = 1
    for i in range(len(sl)):
        slices = sl['number_of_slices'][i]
        cube_name = sl['Cube_id'][i]
        print(sl['Cube_id'][i], sl['number_of_slices'][i])
        cube = np.empty((slices, ax_sz, ax_sz))
        cube_w = np.empty((slices, ax_sz, ax_sz))
        for j in range(slices):
            hdul = fits.open(folder + 'PADDED_%s_image_c%s.%04d.resamp.fits'%(field,chip,idex))
            hdul_w = fits.open(folder + 'PADDED_%s_image_c%s.%04d.resamp.weight.fits'%(field,chip,idex))
            print('PADDED_%s_image_c%s.%04d.resamp.fits'%(field,chip,idex))
            cube[j,:,:] = hdul[0].data    
            cube_w[j,:,:] = hdul_w[0].data    
            idex +=1
        fits.writeto(gd_folder + 'cube%s_gd.fits'%(cube_name),cube, overwrite=True)
        fits.writeto(gd_folder + 'cube%s_gd_w.fits'%(cube_name),cube_w, overwrite=True)
        
# %%
# Erases all the intermediate products.
# for file in glob.glob(os.path.join(folder, "PADDED*.fits")) + glob.glob(os.path.join(folder, "70*.fits")):
#     os.remove(file)
#     print(f"Deleted file: {file}")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        