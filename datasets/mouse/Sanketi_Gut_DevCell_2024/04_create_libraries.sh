#!/bin/bash

project_id='sanketi_2024'
workingdir="$HOME/scratch/ngs/${project_id}"
libraries_dir="${workingdir}/scripts/libraries"

# Check if the directory exists, if yes, remove it
if [ -d "$libraries_dir" ]; then
    rm -rf "$libraries_dir"
fi

# Create the directory again to ensure it's empty
mkdir -p "$libraries_dir"

for folder in ${workingdir}/fastq/*; do
    folder_name=$(basename "$folder")

    if [[ $folder_name == "logs" ]]; then
        continue
    fi

    while IFS=, read -r SRR assay stage tissue strain; do
        if [[ "${folder_name}" == "${SRR}" ]]; then
            break
        fi
    done < "${workingdir}/scripts/sanketi_2024_metadata.csv"

    library_name="${assay}_${strain}_${tissue}_${stage}_GEX"
    library_csv_name=$(echo "$library_name" | sed 's/\.5_GEX/_5_GEX/g')
    output_file="${libraries_dir}/${library_csv_name}.csv"

    if [[ ! -e "$output_file" ]]; then
        {
            echo "[gene-expression]"
            echo "reference,/fast/work/groups/ag_romagnani/ref/mm/GRCm38-hardmasked-optimised-arc"
            echo "create-bam,true"
            echo ""
            echo "[libraries]"
            echo "fastq_id,fastqs,feature_types"
         } > "$output_file"
    fi

    echo "${library_name},${workingdir}/fastq/${SRR},Gene Expression" >> "$output_file"

done
