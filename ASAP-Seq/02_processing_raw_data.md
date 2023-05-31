# Processing raw ASAP-seq from the command line
This guide will be tailored towards the Romagnani Lab, who use the BIH-CUBI cluster. If you ever encounter issues, their is a slightly outdated wiki [here](https://bihealth.github.io/bih-cluster/), and a talk page [here](https://hpc-talk.cubi.bihealth.org/). I will use my username `knighto_c` for this tutorial, but you may see it written as `$USER`, which only applies for sessions where you are logged in. You cannot use $USER in a `slurm`-submitted script, you must explicity set your own username.

## Downloading raw data
Once you receive the notification that your samples have been sequenced, the first step is to download it onto the cluster. The `scratch` directory is used for processing large sets of genetic data, and you can find yours here: `/fast/scratch/users/$USER`.

Depending on where your data was sequenced, the process will be different for downloading. Please ask for login credentials as they will not be written here and will be given out personally.

## Running bcl2fastq2
Once you have downloaded your data, I recommend setting up your project directory as follows:  
1. Rename your downloaded folder, with the RunInfo.xml file present using `mv original_name project_id_bcl`  
2. Create a folder for processing scripts with `mkdir -p project_id_scripts`. (the -p 'flag' is to create it only if it isn't already present)  

Following this, you now need to create two files in your `project_id_scripts` folder
### FASTQ shell script
To submit commands to the job-scheduling program `slurm', you need to create a shell file, which you can do with `nano project_id_scripts/fastq.sh`. Here is a template for generating FASTQ from BCL files for ATAC-seq:

```
#!/bin/bash

#SBATCH -o /fast/home/users/knighto_c/work/slurm/%j.err
#SBATCH -e/fast/home/users/knighto_c/work/slurm/%j.out
#SBATCH --ntasks 32
#SBATCH --mem 64000
#SBATCH --time 6:00:00

#SBATCH -J fqatac
#SBATCH -D /fast/home/users/knighto_c/scratch/ngs

export PATH=/fast/home/users/knighto_c/group/work/bin/cellranger-atac-2.1.0/bin:$PATH
source /fast/home/users/knighto_c/work/bin/miniconda3/etc/profile.d/conda.sh
conda activate sctools

project_id=S3816
project_dir=/fast/home/users/knighto_c/scratch/ngs/${project_id}/

cd $project_dir

cellranger-atac mkfastq --id ${project_id}_fastq --run ${project_id}_bcl --csv ${project_dir}/${project_id}_scripts/indices.csv --use-bases-mask Y88n*,I8n*,Y16n*,Y88n*

rm -r ${project_id}_fastq/_* ${project_id}_fastq/MAKE_FASTQS_CS
```


## Running cellranger-atac count

## Running AMULET and cellsnp-lite/vireoSNP

## Running ASAP-to-KITE

