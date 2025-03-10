#!/bin/bash

# Define the base directory
BASE_DIR="/home/data/alvaro/gns_gd/gns2/F20/cubes_aligned/slices"

# Loop over the four chip directories and run in parallel
for chip in chip1 chip2 chip3 chip4; do
    ls "$BASE_DIR/$chip/"*.[0-9][0-9][0-9][0-9].fits | awk -v chip="$chip" '{print "./sextractor.sh", $1, chip}' | sh &
done

# Wait for all background processes to finish
wait

echo "All processes completed!"
