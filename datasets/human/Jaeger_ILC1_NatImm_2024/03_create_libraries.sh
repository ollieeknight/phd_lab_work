#!/bin/bash

project_id='jaeger_2024'
workingdir="$HOME/scratch/ngs/${project_id}"
libraries_dir="${workingdir}/scripts/libraries"

# Check if the directory exists, if yes, remove it
if [ -d "$libraries_dir" ]; then
    rm -rf "$libraries_dir"
fi

# Create the directory again to ensure it's empty
mkdir -p "$libraries_dir"

# Read metadata CSV file line by line
while IFS=, read -r SRR GSM tissue modality n_donors library_number; do
    echo "Processing metadata entry: $SRR"

    found_match=false

    for folder in "${workingdir}/fastq/"*; do
        folder_name=$(basename "$folder")
        if [[ "${folder_name}" == "${SRR}" ]]; then
            found_match=true
            break
        fi
    done

    if [[ "$found_match" == false ]]; then
        echo "No matching folder found for SRR: $SRR"
        continue
    fi

    library_name="GEX_Jaeger_2024_lib${library_number}_${tissue}"
    output_file="${libraries_dir}/${library_name}.csv"
    echo "Creating output file: $output_file"

    if [[ "${modality}" == "GEX" ]]; then

        if [[ ! -f "$output_file" ]]; then
            {
                echo "[gene-expression]"
                echo "reference,/fast/work/groups/ag_romagnani/ref/hs/GRCh38-hardmasked-optimised-arc"
                echo "create-bam,false"
                echo ""
                echo "[feature]"
                echo "reference,${workingdir}/scripts/jaeger_2024_ADT.csv"
                echo ""
                echo "[libraries]"
                echo "fastq_id,fastqs,feature_types"
            } > "$output_file"
            echo "Header added to: $output_file"
        fi

        echo "${library_name}_${modality},${workingdir}/fastq/${SRR},Gene Expression" >> "$output_file"
        echo "Added GEX line to: $output_file"
    elif [[ "${modality}" == "ADT" ]]; then
        echo "${library_name}_${modality},${workingdir}/fastq/${SRR},Antibody Capture" >> "$output_file"
        echo "Added ADT line to: $output_file"
    fi

done < "${workingdir}/scripts/jaeger_2024_metadata.csv"
