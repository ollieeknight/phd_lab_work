#!/bin/bash

project_id='sanketi_2024'

workingdir="$HOME/scratch/ngs/${project_id}"

cd "${workingdir}/fastq" || exit

for folder in ${workingdir}/fastq/*; do
    folder_name=$(basename $folder)

    if [[ $folder_name == "logs" ]]; then
        continue
    fi

    SRR=""
    while IFS=, read -r SRR assay stage tissue strain; do
        if [[ "${folder_name}" == "${SRR}" ]]; then
            break
        fi
    done < "${workingdir}/scripts/sanketi_2024_metadata.csv"

    if [[ -z "${SRR}" ]]; then
        echo "ERROR: No match found for ${folder_name} in metadata"
        continue
    fi

    echo "Renaming ${SRR} to ${assay}_${strain}_${tissue}_${stage}_GEX"
    mv "${SRR}/${SRR}_1.fastq.gz" "${SRR}/${assay}_${strain}_${tissue}_${stage}_GEX_S1_L001_I1_001.fastq.gz"
    mv "${SRR}/${SRR}_2.fastq.gz" "${SRR}/${assay}_${strain}_${tissue}_${stage}_GEX_S1_L001_R1_001.fastq.gz"
    mv "${SRR}/${SRR}_3.fastq.gz" "${SRR}/${assay}_${strain}_${tissue}_${stage}_GEX_S1_L001_R2_001.fastq.gz"

done
