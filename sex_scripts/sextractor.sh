#!/bin/csh
# to run this: ls /path-to-chip-slices/chip1/*.[0-9][0-9][0-9][0-9].fits | awk '{print"./sextractor.sh",$1,"chip1"}' | sh
if ($#argv < 2) then
  echo "Error: Missing required argument for 'chip'."
  echo "Usage: ./sextractor.sh <file> <chip>"
  exit 1
endif

# Assign variables from command-line arguments
set file = $1
set chip = $2

# Extract filenames
set rfile = $file:t
set cat = $rfile:r".cat"
set mask = $rfile:r".weight.fits"
set psf = $rfile:r".psf"
set head = $rfile:r".head" 
set psfxml = $rfile:r"_psfex.xml"
set prepsfexml = $rfile:r"_prepsfex.xml"
set sexxml = $rfile:r"_sex.xml"
   
  set slices = "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/$chip/"
  set out_put =  "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/sex_slices/$chip/"
  set psf_ima = "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/psfex_ima/$chip/"
  set psf_plots = "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/psfex_plots/$chip/"

# run sextractor - first pass - to feed PSFex
    source-extractor ${slices}$rfile -c prepsfex.sex -CATALOG_NAME ${out_put}$cat -SATUR_LEVEL 50000 -XML_NAME ${out_put}$prepsfexml -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ${slices}$mask 

# Compute the PSF with PSFex

    psfex ${out_put}$cat -c default.psfex -XML_NAME ${out_put}$psfxml -CHECKIMAGE_TYPE  CHI,PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS  -CHECKIMAGE_NAME ${psf_ima}chi.fits,${psf_ima}proto.fits,${psf_ima}samp.fits,${psf_ima}resi.fits,${psf_ima}snap.fits      -CHECKPLOT_NAME      ${psf_plots}selfwhm,${psf_plots}fwhm,${psf_plots}ellipticity,${psf_plots}counts,${psf_plots}countfrac,${psf_plots}chi2,${psf_plots}resi


# Extract final photometry

    source-extractor ${slices}$rfile -c defaultpsf.sex -CATALOG_NAME ${out_put}$cat -PSF_NAME ${out_put}$psf -SATUR_LEVEL 50000 -XML_NAME ${out_put}$sexxml -WEIGHT_TYPE MAP_WEIGHT -WEIGHT_IMAGE ${slices}$mask 
    
end

