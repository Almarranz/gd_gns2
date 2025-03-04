#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:43:44 2025

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

# Makes the cubes bak from the aligned (with VVV slices)
field = 20
slices_folder = '/home/data/alvaro/gns_test/F%s/cubes_aligned/slices/'%(field)
sl = Table.read(slices_folder + '%s_cubes_and_slices.txt'%(field), format = 'ascii')

ax_sz = 2048

# chip = 1
for chip in range(3,5):
   
   
    folder = '/home/data/alvaro/gns_test/F%s/cubes_aligned/test_sex/'%(field)
    
    # for n,fi in enumerate(fits_files):
    idex = 1
    for i in range(len(sl)):
        slices = sl['number_of_slices'][i]
        cube_name = sl['Cube_id'][i]
        print(sl['Cube_id'][i], sl['number_of_slices'][i])
        cube = np.empty((slices, ax_sz, ax_sz))
        cube_w = np.empty((slices, ax_sz, ax_sz))
        for j in range(slices):
            hdul = fits.open(slices_folder+ '%s_image_c%s.%04d.fits'%(field,chip,idex))
            hdul_w = fits.open(slices_folder + '%s_image_c%s.%04d.weight.fits'%(field,chip,idex))
            # Create an ImageHDU for each slice with its header
            slice_hdu = fits.ImageHDU(data=hdul[0].data, header=hdul[0].header)
            slice_hdu_weight = fits.ImageHDU(data=hdul_w[0].data, header=hdul_w[0].header)

            print('%s_image_c%s.%04d.resamp.fits'%(field,chip,idex))
            cube[j,:,:] = hdul[0].data    
            cube_w[j,:,:] = hdul_w[0].data    
            idex +=1
        fits.writeto(folder + 'cube%s_aligned.fits'%(cube_name),cube, overwrite=True)
        fits.writeto(folder + 'cube%s_aligned_w.fits'%(cube_name),cube_w, overwrite=True)
        
# %%
from astropy.io import fits
from astropy.table import Table
import numpy as np

field = 20
slices_folder = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/slices/'%(field)
sl = Table.read(slices_folder + '%s_cubes_and_slices.txt' % (field), format='ascii')

ax_sz = 2048

# Iterate over the chips
for chip in range(1, 2):
    folder = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/test_sex/'%(field)
    idex = 1

    # Loop over each cube defined in the text file
    # for i in range(len(sl)):
    for i in range(2):
        slices = sl['number_of_slices'][i]
        cube_name = sl['Cube_id'][i]
        print(sl['Cube_id'][i], sl['number_of_slices'][i])

        # Create an empty HDUList to store the extensions
        hdul_cube = fits.HDUList([fits.PrimaryHDU()])  # Start with a primary HDU
        hdul_weight = fits.HDUList([fits.PrimaryHDU()])  # For the weight cube

        # Loop over slices to add them as individual extensions
        for j in range(slices):
            # Open the image and weight FITS files
            hdul = fits.open(slices_folder + '%s_image_c%s.%04d.fits' % (field, chip, idex))
            hdul_w = fits.open(slices_folder + '%s_image_c%s.%04d.weight.fits' % (field, chip, idex))

            print('%s_image_c%s.%04d.fits' % (field, chip, idex))
            # Add "FILTER = H" to the headers
            hdul[0].header['FILTER'] = 'H'
            hdul_w[0].header['FILTER'] = 'H'

            # Create an ImageHDU for each slice with its header
            slice_hdu = fits.ImageHDU(data=hdul[0].data, header=hdul[0].header)
            slice_hdu_weight = fits.ImageHDU(data=hdul_w[0].data, header=hdul_w[0].header)

            # Append the HDUs to the HDUList
            hdul_cube.append(slice_hdu)
            hdul_weight.append(slice_hdu_weight)

            idex += 1
            hdul_cube.writeto(folder + 'field%s_chip%s_cube_sl_%s.%02d.fits' % (field,chip,j+1,cube_name), overwrite=True)
            hdul_weight.writeto(folder + 'field%s_chip%s_cube_sl_%s.%02d.weight.fits' % (field,chip,j+1,cube_name), overwrite=True)
        # Save the final HDUList as the cube file
        # hdul_cube.writeto(folder + 'cube%s_aligned.fits' % (cube_name), overwrite=True)
        # hdul_weight.writeto(folder + 'cube%s_aligned.weight.fits' % (cube_name), overwrite=True)
        
        hdul_cube.writeto(folder + 'field%s_chip%s_cube.%02d.fits' % (field,chip,cube_name), overwrite=True)
        hdul_weight.writeto(folder + 'field%s_chip%s_cube.%02d.weight.fits' % (field,chip,cube_name), overwrite=True)
