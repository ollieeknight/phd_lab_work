#!/bin/bash

# Get the bed files directly from Twist online

wget https://www.twistbioscience.com/sites/default/files/resources/2019-09/Twist_Exome_RefSeq_targets_hg38.bed
wget https://www.twistbioscience.com/sites/default/files/resources/2020-04/Twist_MitoPanel_chrM_all_hg38_target.bed

# Merge together and sort

cat Twist* > Twist_Exome_RefSeq_MitoPanel_hg32.bed
bedtools sort Twist_Exome_RefSeq_MitoPanel_hg32.bed > Twist_Exome_RefSeq_MitoPanel_hg32_sorted.bed
