#!/bin/bash

# Function to check if bedtools is installed
check_bedtools() {
    if ! command -v bedtools &> /dev/null; then
        echo "Error: 'bedtools' is not installed or not in the PATH"
        exit 1
    fi
}

check_bedtools

# Get the bed files directly from Twist online
wget https://www.twistbioscience.com/sites/default/files/resources/2019-09/Twist_Exome_RefSeq_targets_hg38.bed
wget https://www.twistbioscience.com/sites/default/files/resources/2020-04/Twist_MitoPanel_chrM_all_hg38_target.bed

# Merge together and sort
cat Twist* > Twist_Exome_RefSeq_MitoPanel_hg32.bed
bedtools sort -i Twist_Exome_RefSeq_MitoPanel_hg32.bed > Twist_Exome_RefSeq_MitoPanel_hg32_sorted.bed