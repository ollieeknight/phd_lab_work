#!/bin/bash

project_id='ty_2023'
workingdir="$HOME/scratch/ngs/${project_id}"
libraries_dir="${workingdir}/scripts/libraries"

# Check if the directory exists, if yes, remove it
if [ -d "$libraries_dir" ]; then
    rm -rf "$libraries_dir"
fi

# Create the directory again to ensure it's empty
mkdir -p "$libraries_dir"

# Read metadata CSV file line by line
while IFS=, read -r SRR origin GSM age sex modality timepoint; do
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

    library_name="CITE_Ty_2023_${sex}${age}_${origin}_${timepoint}"
    output_file="${libraries_dir}/${library_name}.csv"
    echo "Creating output file: $output_file"

    if [[ "${modality}" == "GEX" ]]; then

        if [[ ! -f "$output_file" ]]; then
            {
                echo "[gene-expression]"
                echo "reference,/data/cephfs-2/unmirrored/groups/romagnani/work/ref/hs/GRCh38-hardmasked-optimised-arc"
                echo "create-bam,false"
                echo ""
                echo "[feature]"
                echo "reference,${workingdir}/scripts/${project_id}_ADT_list.csv"
                echo ""
                echo "[libraries]"
                echo "fastq_id,fastqs,feature_types"
            } > "$output_file"
            echo "Header added to: $output_file"
        fi

        echo "CITE_Ty_2023_${sex}${age}_${origin}_${timepoint}_${modality},${workingdir}/fastq/${SRR},Gene Expression" >> "$output_file"
        echo "Added GEX line to: $output_file"

    elif [[ "${modality}" == "ADT" ]]; then
        echo "CITE_Ty_2023_${sex}${age}_${origin}_${timepoint}_${modality},${workingdir}/fastq/${SRR},Antibody Capture" >> "$output_file"
        echo "Added ADT line to: $output_file"
    fi

done < "${workingdir}/scripts/${project_id}_metadata.csv"
