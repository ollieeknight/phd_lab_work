# Processing raw ASAP-seq from the command line
This guide will be tailored towards the Romagnani Lab, who use the BIH-CUBI cluster. If you ever encounter issues, their is a slightly outdated wiki [here](https://bihealth.github.io/bih-cluster/), and a talk page [here](https://hpc-talk.cubi.bihealth.org/). I will use my username `knighto_c` for this tutorial, but you may see it written as `$USER`, which only applies for sessions where you are logged in. You cannot use $USER in a `slurm`-submitted script, you must explicity set your own username.

## Downloading raw data
Once you receive the notification that your samples have been sequenced, the first step is to download it onto the cluster. The `scratch` directory is used for processing large sets of genetic data, and you can find yours here: `/fast/scratch/users/$USER`.

Depending on where your data was sequenced, the process will be different for downloading. Please ask for login credentials as they will not be written here and will be given out personally.

## Converting `BCL` files to `FASTQ`
Once you have downloaded your data, I recommend setting up your project directory as follows:  
1. Rename your downloaded folder, with the RunInfo.xml file present using `mv original_name project_id_bcl`  
2. Create a folder for processing scripts with `mkdir -p project_id_scripts`. (the -p 'flag' is to create it only if it isn't already present)  

Following this, you now need to create two files in your `project_id_scripts` folder

### Sample index sheet
To begin, we need to create an input file which `bcl2fast2` will use to to demultiplex each library using the indices we used during library preparation.  
In the below example, the first two lines correspond to two 10x ATAC v1 libraries (with ASAP-seq modifications), one pool split across two GEM wells. The index name here actually corresponds to four indices which helps capture complexity. You can find the full list [here](https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/sequencing/sample-index-sets-for-single-cell-atac).  
The next four lines correspond to four primers, two of which belong to the antibody capture, and two of which belong to hashtag. We will demultiplex them into one FASTQ file for each sample, but you can also demultiplex them into an ADT and HTO file separately to compare if either library needs more sequencing saturation. You can find the indices for these [here](https://www.biolegend.com/en-us/protocols/totalseq-a-antibodies-and-cell-hashing-with-10x-single-cell-3-reagent-kit-v3-3-1-protocol).

```
nano project_id_scripts/indices.csv
```
```shell
lane,sample,index
*,NK_CMVpos_exp3_libA_atac,SI-NA-A11
*,NK_CMVpos_exp3_libB_atac,SI-NA-A12
*,NK_CMVpos_exp3_libA_adt,AATGAGCG
*,NK_CMVpos_exp3_libB_adt,GGAATCTC
*,NK_CMVpos_exp3_libA_adt,CGCTCATT
*,NK_CMVpos_exp3_libB_adt,GAGATTCC
```

### FASTQ shell script

To submit commands to the job-scheduling program `slurm`, you need to create a shell file. Here is a template for generating FASTQ from BCL files for ATAC-seq:
```
nano project_id_scripts/fastq.sh
```
```shell
#!/bin/bash

#SBATCH -o /fast/home/users/knighto_c/work/slurm/%j.err
#SBATCH -e /fast/home/users/knighto_c/work/slurm/%j.out
#SBATCH --ntasks 32
#SBATCH --mem 64000
#SBATCH --time 6:00:00
#SBATCH -J fq_atac
#SBATCH -D /fast/home/users/knighto_c/scratch/ngs

username=knighto_c
project_id=S1234
bases_mask='Y100n*,I8n*,Y16n*,Y100n*'

ln -s /fast/home/groups/ag_romagnani ~/group
export PATH=/fast/home/users/$username/group/work/bin/cellranger-atac-2.1.0/bin:$PATH
source /fast/home/users/$username/work/bin/miniconda3/etc/profile.d/conda.sh
project_dir=/fast/home/users/$username/scratch/ngs/${project_id}/

conda activate bcl_to_fastq

exec > /fast/home/users/$username/scratch/ngs/${project_id}/fastq_creation.log
{
  conda list
  cd $project_dir
  cellranger-atac mkfastq --id ${project_id}_fastq --run ${project_id}_bcl --csv ${project_dir}/${project_id}_scripts/indices.csv --use-bases-mask $bases_mask
  rm -r ${project_id}_fastq/_* ${project_id}_fastq/MAKE_FASTQS_CS
} 2>&1

mv /fast/home/users/$username/scratch/ngs/${project_id}/fastq_creation.log /fast/home/users/$username/scratch/ngs/${project_id}/${project_id}_fastq
```

Let's break down what's going on here.
1. `#SBATCH` lines tell the system load manager, `slurm`, how many CPUs, RAM, and time you want to use a command for. You then define a name for the job and set your working directory.
2. You then set a few variables: your username, the project_id of the sequencing run, and the bases mask for the experiment, which you can find in the 10x sequencing requirements - but the above is correct for ASAP-seq.
3. Then, if one isn't already made, a symbolic link is made in your home directory for the group folder. The directory for the `cellranger-atac` executable is defined, as well as that for `miniconda`, and a variable is set for your working directory.
4. The `conda` environment generated [here](https://github.com/ollieeknight/single_cell_analysis/blob/main/work-environment/conda_environments.md#processing-raw-bcl-files) is activated.
5. `exec` then opens a connection to write a log at a destination, where the versions of programs used are listed. `cellranger-atac mkfastq` then uses the variables you set to generate fastq files, before removing intermediate files and moving that log into the fastq folder.

Once you have saved these files, you can submit the command with `sbatch project_id_scripts/fastq.sh`

## Generating count matrices with `cellranger-atac count`

```shell
#!/bin/bash

#SBATCH -o /fast/home/users/knighto_c/work/slurm/%j.err
#SBATCH -e /fast/home/users/knighto_c/work/slurm/%j.out
#SBATCH --ntasks 64
#SBATCH --mem 128000
#SBATCH --time 24:00:00

#SBATCH -J atac
#SBATCH -D /fast/home/users/knighto_c/scratch/ngs

username=knighto_c
project_id=S1234
ref=/fast/home/users/$username/group/ref/hs/refdata-arc-hardmasked-GRCh38-2.0.1

export PATH=/fast/home/users/$username/group/work/bin/cellranger-atac-2.1.0/bin:$PATH
project_dir=/fast/home/users/$username/scratch/ngs/$project_id
flowcell_id=$(grep -m 1 'Flowcell=' $project_dir/${project_id}_bcl/RTA3.cfg | sed 's/.*Flowcell=//')
project_fastqs=$project_dir/${project_id}_fastq/outs/fastq_path/$flowcell_id

mkdir -p $project_dir/${project_id}_outs && cd $project_dir/${project_id}_outs

sample_ids=(NK_CMVpos_exp3_libA NK_CMVpos_exp3_libB)

for i in "${!sample_ids[@]}"
do
  sample=${sample_ids[$i]}
	
  cellranger-atac count --id $sample --reference $ref --fastqs $project_fastqs --sample ${sample}_atac
	
  mkdir -p $sample/logs
	mv $sample/*.tgz $sample/logs
	rm -r $sample/SC_ATAC_COUNTER_CS/ $sample/_*

done
```

## Detecting doublets and genotypign cells with `AMULET` and `cellsnp-lite`/`vireoSNP`

## Genotyping mitochondrial DNA with `mgatk`

## Counting antibody capture with `ASAP-to-KITE`


