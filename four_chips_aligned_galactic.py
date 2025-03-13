
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
import time

# %%
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import celestial_frame_to_wcs
from astropy.coordinates import SkyCoord
from astropy import units as u
from reproject import reproject_interp
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import warnings
from astropy.io.fits.verify import VerifyWarning

warnings.simplefilter('ignore', VerifyWarning)

pruebas = '/Users/amartinez/Desktop/for_people/for_Herve/pruebas/'
field = 20
folder = '/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/'

all_chips = '/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/all_chips_gal/'
n_fits = glob.glob(folder + 'chip1/20_image_c1*.weight.fits')


total_t= time.time()
# Read images and WCS headers

for n in range(1, len(n_fits)+1):
    
    images = []
    masks = []
    wcs_list = []
    tick = time.time()
# for n in range(1, 2):
    fits_files = [folder + 'chip1/20_image_c1.%04d.fits'%(n), 
                  folder + 'chip2/20_image_c2.%04d.fits'%(n), 
                  folder + 'chip3/20_image_c3.%04d.fits'%(n), 
                  folder + 'chip4/20_image_c4.%04d.fits'%(n)]
    fits_masks = [folder + 'chip1/20_image_c1.%04d.weight.fits'%(n), 
                  folder + 'chip2/20_image_c2.%04d.weight.fits'%(n), 
                  folder + 'chip3/20_image_c3.%04d.weight.fits'%(n), 
                  folder + 'chip4/20_image_c4.%04d.weight.fits'%(n)]
    for file in fits_files:
        
        with fits.open(file) as hdul:
            images.append(hdul[0].data)  # Read image data
            wcs_list.append(WCS(hdul[0].header))  # Read WCS info
            
    for mask in fits_masks:
        
        with fits.open(mask) as hdul:
            masks.append(hdul[0].data)  # Read image data
            
    
    # Load the original header from `c1`
    with fits.open(fits_files[0]) as hdul:
        original_header = hdul[0].header.copy()  # Make a copy of the header
    
    wcs_keys = [
       "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
       "CDELT1", "CDELT2", "CD1_1", "CD1_2", "CD2_1", "CD2_2",
       "LONPOLE", "LATPOLE", "RADESYS", "EQUINOX"]
    
    for key in wcs_keys:
        original_header.remove(key, ignore_missing=True)

    
    
    # Determine the best common WCS in celestial coordinates
    t_wcs = time.time()
    wcs_out, shape_out = find_optimal_celestial_wcs([(img, wcs) for img, wcs in zip(images, wcs_list)])
    t_wcs_2 = time.time()
    # Convert reference point to Galactic coordinates
    center_icrs = SkyCoord(wcs_out.wcs.crval[0], wcs_out.wcs.crval[1], unit="deg", frame="icrs")
    center_galactic = center_icrs.galactic
    
    # Manually define a Galactic WCS
    wcs_galactic = WCS(naxis=2)
    wcs_galactic.wcs.ctype = ['GLON-TAN', 'GLAT-TAN']  # Set coordinate type to Galactic
    wcs_galactic.wcs.crval = [center_galactic.l.deg, center_galactic.b.deg]  # Set reference point in Galactic coords
    wcs_galactic.wcs.crpix = wcs_out.wcs.crpix  # Keep pixel reference
    wcs_galactic.wcs.cdelt = wcs_out.wcs.cdelt  # Maintain pixel scale
    wcs_galactic.wcs.pc = wcs_out.wcs.pc  # Maintain rotation matrix
    
    # Reproject images to the Galactic WCS
    t_re = time.time()
    mosaic, footprint = reproject_and_coadd([(img, wcs) for img, wcs in zip(images, wcs_list)],
                                            wcs_galactic, shape_out=shape_out, reproject_function=reproject_interp)
    t_re_2 = time.time()
    mosaic_m, footprint_m = reproject_and_coadd([(img, wcs) for img, wcs in zip(masks, wcs_list)],
                                            wcs_galactic, shape_out=shape_out, reproject_function=reproject_interp)
    
    # ---- TRIM THE MOSAIC ----
    # Find non-zero pixel bounding box
    nonzero_pixels = np.where(mosaic > 0)
    y_min, y_max = nonzero_pixels[0].min(), nonzero_pixels[0].max()
    x_min, x_max = nonzero_pixels[1].min(), nonzero_pixels[1].max()
    
    # Trim the mosaic and update WCS
    trimmed_mosaic = mosaic[y_min:y_max, x_min:x_max]
    trimmed_mosaic_m = mosaic_m[y_min:y_max, x_min:x_max]
    trimmed_wcs = wcs_galactic.slice((slice(y_min, y_max), slice(x_min, x_max)))
    
    
    final_header = trimmed_wcs.to_header()
    for key, value in original_header.items():
        final_header[key] = value  # Copy non-WCS keywords
    
    # Save the trimmed Galactic-aligned mosaic
    # hdu = fits.PrimaryHDU(data=trimmed_mosaic, header=trimmed_wcs.to_header())
    hdu = fits.PrimaryHDU(data=trimmed_mosaic, header=final_header)
    hdu_m = fits.PrimaryHDU(data=trimmed_mosaic_m, header=final_header)
    
    # hdu.writeto( all_chips+ "%s_image.%04d.fits"%(n), overwrite=True)
    hdu.writeto(all_chips+ "%s_image.%04d.fits"%(field, n), overwrite=True)
    hdu_m.writeto(all_chips+ "%s_image.%04d.weight.fits"%(field, n), overwrite=True)
    
    print("Trimmed mosaic saved as 'mosaic_galactic_trimmed.fits'")
    
    
    # # Save the Galactic-aligned mosaic
    # hdu = fits.PrimaryHDU(data=mosaic, header=wcs_galactic.to_header())
    # hdu.writeto(pruebas + "mosaic_galactic.fits", overwrite=True)
    
    # print("Mosaic aligned with Galactic coordinates saved as 'mosaic_galactic.fits'")
    tock = time.time()
    print('It took %.1f min'%((tock-tick)/60))
    print('WCS took %.1f min'%((t_wcs_2-t_wcs)/60))
    print('Repre took %.1f min'%((t_re_2-t_re)/60))
    

total_t2= time.time()
print('Total time  %.1f h'%((total_t-total_t2)/3600))























