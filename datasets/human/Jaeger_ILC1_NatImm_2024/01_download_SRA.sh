#!/bin/bash

project_id='jaeger_2024'

workingdir="$HOME/scratch/ngs/${project_id}"

mkdir -p "${workingdir}/fastq/logs"

metadata_file="${project_id}_metadata.csv"

# Read each line from the metadata file and process it
while IFS=, read -r SRR GSM tissue modality n_donors library_number; do
    # Skip the header row
    if [[ $SRR == "SRR" ]]; then
        continue
    fi
    echo "Downloading SRA file ${SRR}"
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${SRR}
#SBATCH --output=${workingdir}/fastq/logs/${SRR}.out
#SBATCH --error=${workingdir}/fastq/logs/${SRR}.out
#SBATCH --ntasks=8
#SBATCH --mem=16000
#SBATCH --time=128:00:00

export PATH=${HOME}/work/bin/pigz:\$PATH
export PATH=~/group/work/bin/sratoolkit.3.1.1-centos_linux64/bin:\$PATH

mkdir -p ${workingdir}/fastq/${SRR}/
cd "${workingdir}/fastq/${SRR}"

fasterq-dump ${SRR} -e \$(nproc) -t ${HOME}/scratch/tmp/ --split-files --include-technical
pigz ${workingdir}/fastq/${SRR}/* -p \$(nproc)
EOF
    echo ""
done < "${metadata_file}"
