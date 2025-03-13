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
ims = '/home/data/GNS/2021/H/20/ims/'
# %%
###########IMPORTANT################
# list of raw images and raw images are in /home/data/raw/GNS_2/H/Field/20
# list of cleand cubes is in '/home/data/GNS/2021/H/20/cleaned'
# list of starfinder stars is in '/home/data/GNS/2021/H/20/data'
# mask for all the chip is un '/home/data/GNS/2021/H/20/ims/mask.fits'
####################################


vvv = Table.read(VVV_fol + 'b333-2024-05-13.dat', format = 'ascii')
lista = open(folder + 'list.txt', 'r')
lines = len(lista.readlines())# This is closing lista. I have to reopen it
lista = open(folder + 'list.txt', 'r')
# %%
hdu_m = fits.open(ims + 'mask.fits')
data_m = hdu_m[0].data
dic_sl = {}

image_i = 0
ch_range = [1,2]


with open(cubes_folder + '%s_cubes_and_slices.txt'%(field),'w') as fil:
    fil.write('Cube_id number_of_slices\n')
for li,l in enumerate(lista):
    
    print(30*'*')
    print(f'Working on file {li+1}/{lines}')
    print(30*'*')
    
    orig_header = fits.getheader(folder + l.strip())
    print(li)
    with gzip.open(clean + 'cube%s.fits.gz'%(li+1),'rb') as f_in:
        with open(clean + 'cube%s.fits'%(li+1),'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
            
    # with fits.open(clean + 'cube%s.fits'%(li+1),mode = 'update') as cube:
    #     cube[0].header = orig_header
    #     cube.flush()   
        # cube.writeto(clean + 'cube%s.fits.gz'%(li+1), overwrite = True)
    cube = fits.open(clean + 'cube%s.fits'%(li+1),mode = 'update') 
    cube[0].header = orig_header
    cube.flush()   
    image_data = cube[0].data
    
    if li>0:
        image_i = image_i + dic_sl[f'l{li-1}']
        print(30*'+',f'\nThe index is {image_i}\n',30*'+')
    
    
    
    wcs = WCS(cube[0].header, naxis = 2)
    naxis1 = cube[0].header['NAXIS1']
    naxis2 = cube[0].header['NAXIS2']
    x_off = cube[0].header['CRPIX2']
    y_off = cube[0].header['CRPIX1']
    
    print(f'Cube {li+1} has %s slices'%(cube[0].header['NAXIS3']))
    
    with open(cubes_folder + '%s_cubes_and_slices.txt'%(field),'a') as fil:
        fil.write('%s %s\n'%(li+1,cube[0].header['NAXIS3']))
    
    use_idx = vvv['J']<900
    vvv_data = vvv[use_idx]
    
    vvv_ra = vvv_data['ra']
    vvv_dec = vvv_data['dec']
    # vvv_gal = SkyCoord(ra = vvv_ra, dec = vvv_dec, unit = 'degree').galactic 
    
    xy = wcs.all_world2pix(np.c_[vvv_data['ra'],vvv_data['dec']],1)
    x = xy[:,0]
    y = xy[:,1]
    
    
    vvv_data.add_columns((xy[:,0],xy[:,1]),names = ('x','y'))
    
    mask = ((x >= -x_off) & (x < naxis1-x_off) & (y >= 0) & (y < naxis2))
    
    # Apply the mask to the vvv.txt data
    vvv_overlap = vvv_data[mask]
    
    
    vvv_overlap.sort('J')
    
    x_vvv = vvv_overlap['x']
    y_vvv = vvv_overlap['y']
    # xh = naxis1/2
    
    
    xh = (max(x_vvv) + min(x_vvv))/2
    yh = (max(y_vvv) + min(y_vvv))/2
    #crop list for each chip
    for chip in range(ch_range[0],ch_range[1]):
        cubes_aligned = '/home/data/alvaro/gns_gd/gns2/F%s/cubes_aligned/slices/chip%s/'%(field,chip)

        
        if (chip == 1):
            idx = np.nonzero((x_vvv < xh) & (y_vvv < yh))
        elif (chip == 2):
            idx = np.nonzero((x_vvv > xh) & (y_vvv < yh))
        elif (chip == 3):
            idx = np.nonzero((x_vvv > xh) & (y_vvv > yh))
        elif (chip == 4):
            idx = np.nonzero((x_vvv < xh) & (y_vvv > yh))
        #
        # fig, ax = plt.subplots(1,1)
        # ax.scatter(x_vvv,y_vvv)
        # ax.scatter(x_vvv[idx],y_vvv[idx])
        # ax.axvline(xh, color = 'r')
        
        gns = Table.read(sf_folder + 'dejitter_stars_%s_%s.txt'%(chip,li+1), format = 'ascii')
        #This ara the xy coordinates from starfinder on the original cubes
        x_gns = gns['x']
        y_gns = gns['y']
        
        xy_gns = np.c_[x_gns,y_gns]
        
        xy_vvv = np.vstack((vvv_overlap['x'][idx],vvv_overlap['y'][idx])).T
        # Astroalign does not work when repeated value appears in the same array
        xy_test, idx_u, = np.unique(xy_vvv, axis=0, return_index=True)
        xy_vvv_unique = xy_vvv[np.sort(idx_u)]
        # sys.exit(154)
        # fig, ax = plt.subplots(1,1)
        # ax.scatter(xy_vvv[:,0],xy_vvv[:,1], s =20)
        # ax.scatter(xy_vvv_unique[:,0],xy_vvv_unique[:,1], s =1)
        np.savetxt(pruebas  + f'GNS_vvv_xy_f20_c{chip}.txt',xy_vvv_unique,fmt = '%.8f')
        
        p, (pos_img, pos_img_t) = aa.find_transform(xy_gns, xy_vvv_unique, max_control_points=200)
        print(p)
        sys.exit(168)
        
        # fig, (ax1,ax2) = plt.subplots(1,2)
        # ax1.scatter(pos_img[:,0],pos_img[:,1])
        # ax2.scatter(pos_img_t[:,0],pos_img_t[:,1])
        
        # Find the common VVV stars in the VVV table to retrive the RA and Dec coordinates
        # x_vvv = pos_img_t  # 
        
    
        x_overlap = np.array(vvv_overlap['x'][idx])
        y_overlap = np.array(vvv_overlap['y'][idx])
        
    
        indices = []
        
    
        for coord in pos_img_t:
            x, y = coord
            # Find the indices where both x and y match in vvv_overlap
            index = np.where((x_overlap == x) & (y_overlap == y))[0]
            
            # If a match is found, append the index
            if index.size > 0:
                indices.append(index[0])  # Assuming one-to-one match
            else:
                indices.append(-1)  # Append -1 for unmatched points
        
        # Convert indices list to numpy array for easier handling
        indices = np.array(indices)
        
        # Print or return the indices
        print('Comon stars in the C%s = %s'%(chip,len(indices)))
        vvv_com = vvv_overlap[idx][indices]
        vvv_com.write(pruebas + 'vvv_common_c%s.txt'%(chip), format = 'ascii', overwrite = True)
        
        
        gns_com_xy = Table(pos_img, names = ('x','y'))
        vvv_com_xy = Table(pos_img_t, names = ('x','y'))
        gns_com_xy.write(pruebas + 'gns_c%s_com.txt'%(chip), format = 'ascii',overwrite=True)
        vvv_com_xy.write(pruebas + 'vvv_c%s_com.txt'%(chip), format = 'ascii',overwrite=True)
        
    
        vvv_com_ad = SkyCoord(ra = vvv_com['ra'],dec = vvv_com['dec'],unit = 'degree')
        xy_com = np.vstack((pos_img[:,0],pos_img[:,1]))
        wcs_new = fit_wcs_from_points(xy_com, vvv_com_ad, projection="TAN")
        
        data_cube = cube[0].data
        header = cube[0].header  # Get the header from the corresponding extension
        m_cube = np.stack([data_m]*cube[0].header['NAXIS3'], axis = 0)
        for i in range(cube[0].header['NAXIS3']):
            image_i +=1
            # Get the data and header for each image (extension)
            if chip ==1:
                data = data_cube[i][0:2048,0:2048]# Use i+1 since the first is the PrimaryHDU
                wcs_header = wcs_new.to_header()                
                m_data = m_cube[i][0:2048,0:2048]
                # for card in wcs_header.cards:
                #     # If the key already exists, it will be replaced, otherwise it will be added
                #     header[card.keyword] = card.value
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value

                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.fits'%(field,chip,image_i),
                             data = data, header = header, overwrite= True)
                # fits.writeto(cubes_aligned + '%s_mask_c%s.%04d.fits'%(field,chip,image_i),
                #              data = m_data, header = header, overwrite= True)
                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.weight.fits'%(field,chip,image_i),
                             data = m_data, header = header, overwrite= True)
            if chip == 2:
                data = data_cube[i][0:2048,2048:]
                wcs_header = wcs_new.to_header()
                m_data = m_cube[i][0:2048,2048:]
                
                # for card in wcs_header.cards:
                #     # If the key already exists, it will be replaced, otherwise it will be added
                #     header[card.keyword] = card.value
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value

                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.fits'%(field,chip,image_i),
                             data = data, header = header, overwrite= True)
                # fits.writeto(cubes_aligned + '%s_mask_c%s.%04d.fits'%(field,chip,image_i),
                #              data = m_data, header = header, overwrite= True)
                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.weight.fits'%(field,chip,image_i),
                             data = m_data, header = header, overwrite= True)
           
            if chip == 3:
                data = data_cube[i][2048:,2048:]
                wcs_header = wcs_new.to_header()
                
                m_data = m_cube[i][2048:,2048:]
                # for card in wcs_header.cards:
                #     # If the key already exists, it will be replaced, otherwise it will be added
                #     header[card.keyword] = card.value
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value

                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.fits'%(field,chip,image_i),
                             data = data, header = header, overwrite= True)
                # fits.writeto(cubes_aligned + '%s_mask_c%s.%04d.fits'%(field,chip,image_i),
                #              data = m_data, header = header, overwrite= True)
                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.weight.fits'%(field,chip,image_i),
                             data = m_data, header = header, overwrite= True)
              
               
            if chip == 4:
                data = data_cube[i][2048:,0:2048]
                wcs_header = wcs_new.to_header()
                
                m_data = m_cube[i][2048:,0:2048]
                # for card in wcs_header.cards:
                #     # If the key already exists, it will be replaced, otherwise it will be added
                #     header[card.keyword] = card.value
                for card in wcs_header.cards:
                    # Only replace the existing keyword in the header if it exists
                    if card.keyword in header:
                        header[card.keyword] = card.value

                # Create an ImageHDU with the data and the corresponding header
                # image_hdu = fits.ImageHDU(data=data, header=wcs_header)
                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.fits'%(field,chip,image_i),
                             data = data, header = header, overwrite= True)
                # fits.writeto(cubes_aligned + '%s_mask_c%s.%04d.fits'%(field,chip,image_i),
                #              data = m_data, header = header, overwrite= True)
                fits.writeto(cubes_aligned + '%s_image_c%s.%04d.weight.fits'%(field,chip,image_i),
                             data = m_data, header = header, overwrite= True)
        image_i -=  cube[0].header['NAXIS3']
    print(30*'+',f'\nAfter cube{li}, index ={image_i} \n',30*'+')
    dic_sl['l%s'%(li)] = cube[0].header['NAXIS3']
    os.remove(clean + 'cube%s.fits'%(li+1))
   
    # if li == 0:
    #     break    
sys.exit(296)
# %%
# #MISSFITS

# for chip in range(ch_range[0],ch_range[1]):
#     command = ['missfits', cubes_aligned + '%s_image_c%s'%(field,chip), '-c', 'conf.missfits']
    
#     try:
#         # Run the command
        
#         result = subprocess.run(command, check=True)
#         # Print standard output and error
#         print("Command Output:")
#         print(result.stdout)
#         print("Command Error (if any):")
#         print(result.stderr)
    
#     except subprocess.CalledProcessError as e:
#         # Handle errors
#         print(f"Error: {e}")
#         print(f"Standard Output: {e.stdout}")
#         print(f"Standard Error: {e.stderr}")


# # %%
# for chip in range(ch_range[0],ch_range[1]):
#     command = ['missfits', cubes_aligned + '%s_mask_c%s'%(field,chip), '-c', 'conf.missfits']
    
#     try:
#         # Run the command
        
#         result = subprocess.run(command, check=True)
#         # Print standard output and error
#         print("Command Output:")
#         print(result.stdout)
#         print("Command Error (if any):")
#         print(result.stderr)
    
#     except subprocess.CalledProcessError as e:
#         # Handle errors
#         print(f"Error: {e}")
#         print(f"Standard Output: {e.stdout}")
#         print(f"Standard Error: {e.stderr}")
# answer = input('Did you check the generated cubes? \nIf yes and everthing is alright, type "y".\nOtherwise do so and elimate bad slices (see delete_bad.py')
# if answer == 'y':
#     print('Running sextractor, scamp and SWarp')
# else:
#     sys.exit(6*'😡'+'\nCHECK THEM!!!\n'+6*'😡')








        
