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
# %%

locat = '/Volumes/teabag-alvaro/'
# locat = '/home/data/alvaro/'

field = 20
# chip = 'all_chips'
chip = 'chip1'



vvv_f = locat + 'gns_gd/gns2/VVV_fields/'
s_sli = locat + f'gns_gd/gns2/F{field}/cubes_aligned/slices/sex_slices/{chip}/'




hdul_cat = fits.open(s_sli +'2MASS_1743-2840_r4.cat')


vvv = Table.read(vvv_f +  f'vvv_f{field}.dat', format = 'ascii')



# Assuming 'vvv' is your existing Astropy Table
# Create a new table with the first three columns
vvv_cat = vvv['ra', 'dec', 'Ks']

# Rename the columns
vvv_cat.rename_column('ra', 'X_WORLD')
vvv_cat.rename_column('dec', 'Y_WORLD')
vvv_cat.rename_column('Ks', 'MAG')

# Set the units for each column
vvv_cat['X_WORLD'].unit = u.degree
vvv_cat['Y_WORLD'].unit = u.degree
vvv_cat['MAG'].unit = u.mag
col_add = ['ERRA_WORLD','ERRB_WORLD','PMALPHA_J2000','PMDELTA_J2000',
           'PMALPHAERR_J2000','PMALPHAERR_J2000','PMDELTAERR_J2000','MAGERR', '>f4']
vvv_cat['OBSDATE'] = 2006
for i in col_add:
    vvv_cat[i] = np.nan



# Display the new table
print(vvv_cat)

# vvv_cat.write(s_sli + 'vvv.cat', format ='ascii')
vvv_cat.write(s_sli + 'vvv.cat', format ='fits', overwrite = True)

# sys.exit()
# %%
vvv_f = fits.open(s_sli +'vvv.cat')


for i,k in enumerate(hdul_cat[2].header):
    print(i)
    # Only replace the existing keyword in the header if it exists
    # print(k,hdul_cat[2].header[k])
    if i > 7:
        try:
            if vvv_f[1].header[k] == hdul_cat[2].header[k]:
                print('SAME',vvv_f[1].header[k])
            else:
                vvv_f[1].header[k] = hdul_cat[2].header[k]
        except:
            vvv_f[1].header[k] = hdul_cat[2].header[k]
            print(f'No {k} card')

vvv_f.writeto(s_sli + 'vvv.cat', overwrite= True)
          
# %%
        
new_cat =  locat + 'gns_gd/VVV_cat_scamp/'
vvv_t = fits.open(new_cat +'vvvdr4.refcat')





























