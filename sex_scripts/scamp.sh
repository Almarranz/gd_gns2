#!/bin/csh

foreach file ($*)
  set rfile = $file:t
  set cat = $rfile:r".cat"
  set mask = $rfile:r".weight.fits"
  set psf = $rfile:r".psf"
  set head = $rfile:r".head" 
  set psfxml = $rfile:r"_psfex.xml"
  set prepsfexml = $rfile:r"_prepsfex.xml"
  set sexxml = $rfile:r"_sex.xml"

  # Define the full directory path for each plot type
  
   
  set chip = "all_chips"
  set slices = "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/$chip/"
  set out_put =  "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/sex_slices/$chip/"
  set plot_dir = "/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices/scamp_plots/$chip/"

  # Run SCAMP with the modified plot settings, ensuring full path for each plot
  scamp ${out_put}$cat -c ./scamp.conf -CHECKPLOT_NAME ${plot_dir}fgroups,${plot_dir}distort,${plot_dir}astr_interror2d,${plot_dir}astr_interror1d,${plot_dir}astr_referror2d,${plot_dir}astr_referror1d,${plot_dir}astr_chi2,${plot_dir}psphot_error -CHECKPLOT_DEV PDF -MERGEDOUTCAT_NAME   ${out_put}merged.ocat  -XML_NAME   ${out_put}scamp.xml -REFOUT_CATPATH ${out_put}
end
