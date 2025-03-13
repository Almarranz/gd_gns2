#!/bin/bash

cd /home/data/alvaro/gns_gd/gns2/F11/pruebas/

for file in GC_H_F11.H.2022-06-10T05_15_54_*.*.fits; do 
    base="${file%_*.*.fits}"  # Extracts the base filename
    id=$(echo "$file" | awk -F'[_.]' '{print $(NF-2)}')  # Extracts the identifier
    id_padded=$(printf "%03d" "$id")  # Formats with leading zeros
    mv "$file" "${base}_4.${id_padded}.fits"
done
