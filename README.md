# Geometric Distortions in GNS

## Pipeline Overview

There is not yet a defined pipeline. Below are the names and descriptions of some scripts used in the process.

## Prepare the Data

- This pipeline uses reduced FITS images from the GNS2 pipeline. These images include flat, sky, etc., and are reduced.
- The images are divided by chips and have no astrometric solution.

### Delete Bad Images

- `gns_cubes_for_inspection.py`: Groups all slices from the cleaned cubes into four cubes (one per chip). This approach is necessary for visual inspection. Note the positions of bad slices and pass this list to the next script.
- `delete_bad_from_cubes.py`: Deletes bad slices in the cubes. Feed it with a list of bad slices.

  - **Note**: Consider moving this procedure directly to the GNS2 pipeline so that cubes from the `cleaned` directory are actually cleaned.
  - **Note 2**: You can combine this inspection with `maxitract`. Although it doesn't work perfectly, it can help detect elusive bad slices. See `maxitract_chips.sh` and `maxitract_bad_slices.py`.

    > `maxitract_chips.sh`: Runs `maxitract` on FITS slices for different chips and renames the output file to `maxitract.output`. This script can be found in the GNS1 folder.
    > `maxitract_bad_slices.py`: Combines all `maxitract` output files into a single list. This is convenient for comparing bad slices selected by `maxitract` with those selected by the operator.

### Initial Geometric Solution

- For the Astromatic machinery to work, the images need an astrometric solution with an accuracy of a few arcseconds.
- `cubes_for_gd_corrections.fits`: Provides these files with an astrometric solution. It finds matches between the VVVx catalog and the stars present in each chip (using `astroalign` and `sextractor`).


- Runs `missfits`on cubes to divide them into slices. Then make a .list file containing the paths and names of the slices and feed with it `maxitrack`
- **% missfits -d > default.missfits** generates de default configurataion file. 
- On the configuration file you can set `OUTFILE_TYPE  SPLIT` and `SAVE_TYPE NEW`. This would slice the MEF and conserve the original file, while saving the slices under a new name.
- Make a .list files with then names and tha absolutes paths of the new genareted *miss* slices: ** ls *.miss.fits | xargs realpath > part[1,1]_c[1,2,3,4]_fits.list **


## Astromatic

### Photometry
- `./sextractor.sh`
- We have to run `sextrator` to a first estimate of the psf
- Then we have to run `PSFex` to extract a field variable psf
- Finally we run `sextracto` agaim, feeded with the vairable psf

### SCAMP
- `./scamp.sh`
- We use `scamp` to calculate the variability of the *pixel scale* acrros the chip, the geomtric distortion
- It has to be use in the individuals fits files.



















# WARNING
---
This is the old readme, not useful anymore!!!. 

---
## Steps to Estimate and Apply Geometric Distortions (GD) on GNS2

---

### 0. Prepare Cubes for GD Corrections
- Run `cubes_for_gd_corrections.py` to group all the cleaned cubes (reduced for sky, flat-fielded, and cleaned for cosmic rays) into four FITS files, one for each chip.
- This script automatically runs `missfits`, which groups all the slices into a cube and deletes the individual slices afterward.
- The generated cubes will have the correct format for use with Astromatic software. **Inspect the generated cubes carefully**.

#### 0.1 Optional: Remove Bad Frames
- Use `delete_bad.py` if additional bad slices are detected during inspection. In principle, the cubes used in `cubes_for_gd_corrections.py` should already be cleaned of bad frames.
- The script moves the bad slices to a `bad_slices` folder for each field. (**Send these folders to Hervé**.)

#### 0.2 Optional: Check Bad Frames:
- `cubes_for_checking.py`it combined slices in a cube for ocular inspection. It uses `missfits` software.
- `rename_fits.sh` rename some files that Hervé send me to adecuate them to the format required by `missfits`
---

### Astromatic.py
- This script runs all three steps required for GD corrections (#1, #2, and #3), namely: `sextractor`, `scamp`, and `swarp`.

---

### Alternative: Running Steps Separately
If desired, you can run `sextractor`, `scamp`, and `swarp` separately as described below:

#### 1. Run Source Extractor
- Execute `source-extractor cube.fits -c default_c[1,2,3,4].sex`. This will generate a **MEF (Multi-Extension FITS)** file with different pointings, with an extension for each image.
- Ensure that this cube was created using the `cubes_for_gd_corrections.py` script.
- Edit the `default.param` file to select the variables required by SCAMP. A list of all available SExtractor variables can be found in `default.sex`.

#### 2. Run SCAMP
- Run the command: `scamp sextractor.cat -c scamp.conf`.
- SExtractor will generate a catalog to feed into SCAMP, which calculates the geometric distortion. SCAMP will generate a `header.head` file.

#### 3. Run SWarp
- Run `SWarp cube.fits`. 
- The `.head` file generated by SCAMP is used to modify the original `cube.fits`. Ensure this `.head` file is correctly referenced in the configuration.

---

### 4. Apply Axis Correction
- Run `axis_correction.py` to equalize the axes in the FITS files generated by SWarp. This step updates the WCS information accordingly.

---

Follow these steps to properly apply geometric distortion corrections to GNS2 data.

