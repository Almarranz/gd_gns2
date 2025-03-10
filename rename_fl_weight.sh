#!/bin/bash

# Loop over each file matching the pattern
for file in 20_image_c*.weight.fits; do
  # Extract the file ID and chip ID from the filename
  if [[ $file =~ 20_image_c([0-9]+)\.([0-9]{4})\.weight.fits ]]; then
    chip_number=${BASH_REMATCH[1]}
    file_id=${BASH_REMATCH[2]}

    # Format the chip number as a two-digit number
    formatted_chip_number=$(printf "%02d" $chip_number)

    # Construct the new filename
    new_file="20_image_${file_id}.${formatted_chip_number}.weight.fits"

    # Rename the file
    mv "$file" "$new_file"

    # Output the renaming operation for verification
    echo "Renamed '$file' to '$new_file'"
  fi
done
