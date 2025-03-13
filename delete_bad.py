#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 15:38:03 2024

@author: alvaro
"""
from astropy.io import fits
import os

# =============================================================================
# THE FIRST cell removes bad extension from the cube
# =============================================================================
field = 11
# cubes_aligned = f'/home/data/alvaro/gns_gd/gns2/F{field}/cubes_aligned/slices/'
cubes_aligned = f'/home/data/alvaro/gns_gd/gns2/F{field}/pruebas/cubes/'
bad_folder = f'/home/data/alvaro/gns_gd/gns2/bad_slices/F{field}/'

# Ensure the bad_folder exists
os.makedirs(bad_folder, exist_ok=True)

# Iterate over the chips
for chip in range(1, 5):
    # fits_file = cubes_aligned + f"/chip{chip}/{field}_image_c{chip}.fits"
    fits_file = cubes_aligned + f"GC_H_F11.H.2022-06-10T05_15_54_{chip}.fits"
    # extensions_to_remove = [42,	47,	82,	222,	223,	241,	277,	386,	387]
    extensions_to_remove = [56,	57,	58,	324,	387,	388,	401,	415,	429,	456,]

    # Sort extensions in descending order for safe removal
    extensions_to_remove.sort(reverse=True)

    # Open the FITS file
    with fits.open(fits_file, mode='update') as hdulist:
        for ext in extensions_to_remove:
            # Create a new FITS file for the bad slice
            bad_slice_path = os.path.join(bad_folder, f"bad_slice_{chip}_ext{ext}.fits")
            hdu = hdulist[ext]

            # Save the bad slice to the bad folder
            hdu.writeto(bad_slice_path, overwrite=True)

            # Remove the extension from the original FITS file
            # del hdulist[ext]

        # Save changes to the FITS file
        hdulist.flush()

    print(f"Processed chip {chip}, removed extensions: {extensions_to_remove}")

# from astropy.io import fits
# =============================================================================
# for chip in range(1, 5):
#     fits_file = cubes_aligned + f"/chip{chip}/{field}_image_c{chip}.weight.fits"
#     extensions_to_remove = [42,	47,	82,	222,	223,	241,	277,	386,	387]
# 
#     # Sort extensions in descending order for safe removal
#     extensions_to_remove.sort(reverse=True)
# 
#     # Open the FITS file
#     with fits.open(fits_file, mode='update') as hdulist:
#         for ext in extensions_to_remove:
#             # Create a new FITS file for the bad slice
#             bad_slice_path = os.path.join(bad_folder, f"bad_slice_{chip}_ext{ext}.fits")
#             hdu = hdulist[ext]
# 
#             # Save the bad slice to the bad folder
#             hdu.writeto(bad_slice_path, overwrite=True)
# 
#             # Remove the extension from the original FITS file
#             del hdulist[ext]
# 
#         # Save changes to the FITS file
#         hdulist.flush()
# 
#     print(f"Processed chip {chip}, removed extensions: {extensions_to_remove}")
# =============================================================================
    

#%%

# =============================================================================
# THE SECOND cell remove fits files 
# =============================================================================
# =============================================================================
# import os
# field = 20
# cubes_aligned = f'/home/data/alvaro/gns_gd/gns2/F{field}/cubes_aligned/slices/'
# 
# # List of file indices to delete
# files_to_delete = [42, 47, 82, 222, 223, 241, 277, 386, 387]
# for chip in range(1,5):
#     
# # Loop over the list and delete the files
#     for num in files_to_delete:
#         folder = cubes_aligned + f"chip{chip}/"
#         file1 = os.path.join(folder, f"20_image_c{chip}.{num:04d}.fits")
#         file2 = os.path.join(folder, f"20_image_c{chip}.{num:04d}.weight.fits")
#     
#         for file in [file1,file2]:
#             if os.path.exists(file):
#                 os.remove(file)
#                 print(f"Deleted: {file}")
#             else:
#                 print(f"File not found: {file}")
# =============================================================================
