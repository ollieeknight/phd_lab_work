#!/bin/bash

# Define directories
fastq_dir="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/EAE/fastq"
output_dir="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/EAE/mixcr"
log_dir="/data/cephfs-1/scratch/groups/romagnani/users/knighto_c/ngs/EAE/logs"

# Create output and log directories if they don't exist
mkdir -p "${output_dir}"
mkdir -p "${log_dir}"

# Loop through each R1 file
for R1 in "$fastq_dir"/*_R1_001.fastq.gz; do
    # Extract sample name by removing _R1_001.fastq.gz
    sample=$(basename "$R1" _R1_001.fastq.gz)

    # Remove the _Sxx_Lxxx part from the sample name
    clean_sample=$(echo "$sample" | sed -E 's/_S[0-9]+_L[0-9]+//')

    # Define the corresponding R2 file
    R2="$fastq_dir/${sample}_R2_001.fastq.gz"

    # Create output directory for the clean sample
    mkdir -p "$output_dir/${clean_sample}"

    # Submit the SLURM job directly using sbatch and a heredoc
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=MiXCR
#SBATCH --output=${log_dir}/${clean_sample}.out
#SBATCH --error=${log_dir}/${clean_sample}.out
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=96:00:00

# Activate the Conda environment
source ${HOME}/work/bin/miniforge3/etc/profile.d/conda.sh
conda activate mixcr

# Run MiXCR command for this sample
mixcr analyze rna-seq --species mmu --rna -t 32 "$R1" "$R2" "$output_dir/${clean_sample}/"
EOF

done

