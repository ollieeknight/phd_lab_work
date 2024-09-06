#!/bin/bash

project_id="ty_2023"

workingdir="$HOME/scratch/ngs/${project_id}"

mkdir -p "${workingdir}/fastq/logs"

while IFS=, read -r SRR origin GSM age sex modality timepoint; do
    if [[ "$SRR" == "SRR" ]]; then
        continue
    fi

    echo "Submitting download job for $SRR"
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name ${SRR}
#SBATCH --output "${workingdir}/fastq/logs/${SRR}.out"
#SBATCH --error "${workingdir}/fastq/logs/${SRR}.out"
#SBATCH --ntasks=2
#SBATCH --mem=8000
#SBATCH --time=96:00:00
export PATH=${HOME}/work/bin/pigz:\$PATH
export PATH=~/group/work/bin/sratoolkit.3.1.1-centos_linux64/bin:\$PATH

mkdir -p ${workingdir}/fastq/${SRR}/
cd "${workingdir}/fastq/${SRR}"

fasterq-dump ${SRR} -e \$(nproc) -t ${HOME}/scratch/tmp/ --split-files --include-technical
pigz ${workingdir}/fastq/${SRR}/* -p \$(nproc)

EOF
done < "${workingdir}/scripts/${project_id}_metadata.csv"
