from astropy.io import fits
import os
import sys
import glob



# Path where the FITS files are stored
# local = '/home/data/alvaro/'
local = '/Volumes/teabag-alvaro/'

for chip in range(1,2):
    # path = local  + '/gns_gd/gns2/F20/cubes_aligned/slices/chip%s/'%(chip)
    path = local  + 'HB_red/pruebas/GNS2/F4/raw/tmp/'
   
    # n_file = glob.glob(local  + 'HB_red/pruebas/GNS2/F4/raw/tmp/*.fits')
    n_file = glob.glob('/home/data/alvaro/HB_red/pruebas/GNS2/F4/raw/tmp/*.fits')
    # print(n_file)
   
    # Loop over the series of FITS files
    for file in n_file:  # adjust the range based on how many files you have
    # for i in range(1, len(n_file)+1):  # adjust the range based on how many files you have
        # filename = f"20_image_c{chip}.{i:04d}.weight.fits"  # Assuming 4 digits (0001, 0002, ..., 0020)
        # filepath = os.path.join(path, filename)
        
        
        try:
            # Open the FITS file
            # with fits.open(filepath, mode='update') as hdul:
            with fits.open(file, mode='update') as hdul:
                # Add the 'FILTER' keyword with value 'H' to the header of the primary HDU
                hdr = hdul[0].header
                hdr['FILTER'] = 'H'
                print(f"Added FILTER keyword to {file}")
                
                # Save changes
                hdul.flush()
        
        except Exception as e:
            print(f"Error processing {file}: {e}")
        
    
    # # Loop over the series of FITS files
    # for i in range(1, len(n_file)+1):  # adjust the range based on how many files you have
    #     filename = f"20_image_c{chip}.{i:04d}.fits"  # Assuming 4 digits (0001, 0002, ..., 0020)
    #     filepath = os.path.join(path, filename)
        
    #     try:
    #         # Open the FITS file
    #         with fits.open(filepath, mode='update') as hdul:
    #             # Add the 'FILTER' keyword with value 'H' to the header of the primary HDU
    #             hdr = hdul[0].header
    #             hdr['FILTER'] = 'H'
    #             print(f"Added FILTER keyword to {filename}")
                
    #             # Save changes
    #             hdul.flush()
        
    #     except Exception as e:
            # print(f"Error processing {filename}: {e}")
            
