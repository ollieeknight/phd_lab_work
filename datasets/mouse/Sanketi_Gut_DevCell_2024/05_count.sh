#!/bin/bash

project_id='sanketi_2024'
workingdir="$HOME/scratch/ngs/${project_id}"

# Create necessary directories
mkdir -p "${workingdir}/outs/logs/"
cd "${workingdir}/outs"

# Loop through library CSV files
for library_csv in "${workingdir}/scripts/libraries/"*; do
library_id=$(basename "${library_csv%.*}")

    # Check for GEX libraries
    if [[ "${library_id}" == *_GEX ]]; then

        # Check if the output directory does not exist
        if [[ ! -d "${workingdir}/outs/${library_id}/outs" || ! -d "${workingdir}/outs/${library_id}/SC_MULTI_CS" ]]; then

            # Submit cellranger job
            job_id=$(sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=${library_id}
#SBATCH --output="${workingdir}/outs/logs/${library_id}_cellranger.out"
#SBATCH --error="${workingdir}/outs/logs/${library_id}_cellranger.out"
#SBATCH --ntasks=32
#SBATCH --mem=96000
#SBATCH --time=96:00:00

num_cores=\$(nproc)
container="/fast/scratch/users/knighto_c/tmp/oscar-count_latest.sif"
cd "${workingdir}/outs/"
echo "apptainer run -B /fast,/data "\$container" cellranger multi --id "${library_id}" --csv "${library_csv}" --localcores "\$num_cores""
echo ""
apptainer run -B /fast,/data "\$container" cellranger multi --id "${library_id}" --csv "${library_csv}" --localcores "\$num_cores"
rm -r "${workingdir}/outs/${library_id}/SC_MULTI_CS" "${workingdir}/outs/${library_id}/_"*
EOF
            )

            job_id=$(echo $job_id | awk '{print $4}')

            sbatch --dependency=afterok:$job_id <<EOF
#!/bin/bash
#SBATCH --job-name=${library_id}_cellbender
#SBATCH --output="${workingdir}/outs/logs/${library_id}_cellbender.out"
#SBATCH --error="${workingdir}/outs/logs/${library_id}_cellbender.out"
#SBATCH --ntasks=1
#SBATCH --partition="gpu"
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128000
#SBATCH --time=24:00:00

cd "${workingdir}/outs/${library_id}"
mkdir -p cellbender
container="/fast/scratch/users/knighto_c/tmp/oscar-qc_latest.sif"
apptainer run --nv -B /fast,/data "\$container" cellbender remove-background --cuda --input outs/multi/count/raw_feature_bc_matrix.h5 --output cellbender/output.h5
rm ckpt.tar.gz
EOF
        fi
    fi
done
