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
*,ASAP_NKG2C_exp3_libA_atac,SI-NA-A11
*,ASAP_NKG2C_exp3_libB_atac,SI-NA-A12
*,ASAP_NKG2C_exp3_libA_adt,AATGAGCG
*,ASAP_NKG2C_exp3_libB_adt,GGAATCTC
*,ASAP_NKG2C_exp3_libA_adt,CGCTCATT
*,ASAP_NKG2C_exp3_libB_adt,GAGATTCC
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
#SBATCH -J K0006
#SBATCH -D /fast/home/users/knighto_c/scratch/ngs

project_id=S3816
bases_mask='Y88n*,I8n*,Y16n*,Y88n*'

export PATH=~/group/work/bin/cellranger-atac-2.1.0/bin:$PATH
source ~/work/bin/miniconda3/etc/profile.d/conda.sh
project_dir=~/scratch/ngs/${project_id}/

cd $project_dir

conda activate bcl_to_fastq

exec > ~/scratch/ngs/${project_id}/fastq_creation.log
{
        conda list

        cellranger-atac mkfastq --id ${project_id}_fastq --run ${project_id}_bcl --csv ${project_id}_scripts/indices.csv --use-bases-mask $bases_mask

        rm -r ${project_id}_fastq/_* ${project_id}_fastq/MAKE_FASTQS_CS
        mv $project_dir/fastq_creation.log $project_dir/${project_id}_fastq/fastq_creation.log
} 2>&1
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
#SBATCH --time 72:00:00
#SBATCH -J S3816
#SBATCH -D /fast/home/users/knighto_c/scratch/ngs

project_id=S3816
ref=~/group/work/ref/hs/arc-hardmasked-optimised-GRCh38

export PATH=~/group/work/bin/cellranger-atac-2.1.0/bin:$PATH

project_dir=~/scratch/ngs/${project_id}/

flowcell_id=$(grep -m 1 'Flowcell=' ~/scratch/ngs/$project_id/${project_id}_bcl/RTA3.cfg | sed 's/.*Flowcell=//')
project_fastqs=~/scratch/ngs/$project_id/${project_id}_fastq/outs/fastq_path/$flowcell_id
echo $project_fastqs

cd $project_dir

mkdir -p ${project_id}_outs && cd ${project_id}_outs

sample_ids=(ASAP_NKG2C_exp3_libA ASAP_NKG2C_exp3_libB)

for i in "${!sample_ids[@]}"
do
        sample=${sample_ids[$i]}

        cellranger-atac count --id $sample --reference $ref --fastqs $project_fastqs --sample ${sample}_atac

        mkdir -p $sample/logs
        mv $sample/*.tgz $sample/logs
        rm -r $sample/SC_ATAC_COUNTER_CS/ $sample/_*
done
```

## Post-processing
For the next few steps, we actually run everything in one submitted script. I have separated these by each function, however be sure to paste all of it together if you would like to also perform everything in one step. Critical things to remember here including ensure [all conda environments are created](https://github.com/ollieeknight/single_cell_analysis/blob/main/work-environment/conda_environments.md), and that the naming is conserved between it.  
```shell
#!/bin/bash
#SBATCH -o /fast/home/users/knighto_c/work/slurm/%j.err
#SBATCH -e /fast/home/users/knighto_c/work/slurm/%j.out
#SBATCH --ntasks 32
#SBATCH --mem 128000
#SBATCH --time 24:00:00

#SBATCH -J geno
#SBATCH -D /fast/home/users/knighto_c/scratch/ngs

project_id=S1234

source ~/work/bin/miniconda3/etc/profile.d/conda.sh
snp_ref=~/group/work/ref/vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz

project_dir=~/scratch/ngs/${project_id}

cd ${project_dir}/${project_id}_outs

sample_ids=(NK_CMVpos_exp3_libA NK_CMVpos_exp3_libB)
donors_multiplexed=(2 2 2)

```
### Detecting doublets and genotyping cells with `AMULET` and `cellsnp-lite`/`vireoSNP`
```shell
for i in ${!sample_ids[@]}
do
        sample=${sample_ids[$i]}
        donors=${donors_multiplexed[$i]}

        mkdir -p $sample/logs

        exec > $sample/logs/genotype_donors.log
        {
                conda activate donor_genotyping
                conda list
                if [ -d $sample/genotype_peaks ]; then rm -rf $sample/genotype_peaks; fi
                mkdir -p $sample/genotype_peaks
                cellsnp-lite -s $sample/outs/possorted_bam.bam -b $sample/outs/filtered_peak_bc_matrix/barcodes.tsv -O $sample/genotype_peaks -R $snp_ref -p 20 --minMAF 0.1 --minCOUNT 20 --gzip -p 32 --UMItag None
                vireo -c $sample/genotype_peaks -o $sample/genotype_peaks -N $donors -p 32

                conda deactivate
        } 2>&1

        exec > $sample/logs/genotype_mito.log
        {
                conda activate mitochondrial_genotyping
                conda list

                if [ -d $sample/genotype_mito ]; then rm -rf $sample/genotype_mito; fi
                mkdir -p $sample/genotype_mito
                mgatk tenx -i $sample/outs/possorted_bam.bam -n $sample -o $sample/genotype_mito -c 8 -bt CB -b $sample/outs/filtered_peak_bc_matrix/barcodes.tsv

                rm -r .snakemake

                conda deactivate
        } 2>&1

        exec > $sample/logs/amulet_doublets.log
        {
                conda activate amulet
                conda list

                if [ -d $sample/amulet ]; then rm -rf $sample/amulet; fi
                mkdir -p $sample/amulet
                bash ~/group/work/bin/amulet/AMULET.sh $sample/outs/fragments.tsv.gz $sample/outs/singlecell.csv ~/group/work/bin/amulet/human_autosomes.txt ~/group/work/bin/amulet/RestrictionRepeatLists/restrictionlist_repeats_segdups_rmsk_hg38.bed $sample/amulet/ ~/group/work/bin/amulet/

                conda deactivate
        } 2>&1
done
```

### Counting antibody capture with `ASAP-to-KITE`
```shell
project_id=K0004
threads=32

source ~/work/bin/miniconda3/etc/profile.d/conda.sh
featuremap=~/group/work/bin/asap/featuremap.py
asap_to_kite=~/group/work/bin/asap/asap_to_kite_v2.py
atac_whitelist=~/group/work/bin/whitelists/737K-cratac-v1.txt
project_dir=~/scratch/ngs/$project_id
project_fastqs=$project_dir/${project_id}_fastq/outs/fastq_path
mkdir -p $project_fastqs/asap_fastqs

cd $project_dir/${project_id}_outs

conda activate adt_count

sample_ids=(NK_CMVpos_exp3_libA NK_CMVpos_exp3_libB)

mkdir -p adt_index
python $featuremap $project_dir/${project_id}_scripts/asap_adt.csv --t2g adt_index/FeaturesMismatch.t2g --fa adt_index/FeaturesMismatch.fa --header --quiet
kallisto index -i adt_index/FeaturesMismatch.idx -k 15 adt_index/FeaturesMismatch.fa -o adt_index/

for i in "${!sample_ids[@]}"
do
        sample=${sample_ids[$i]}
        mkdir -p "$sample/logs"

        {
                exec > "$sample/logs/adt_counting.log"

                conda list

                python "$asap_to_kite" -f "$project_fastqs" -s "${sample}_adt*" -o "$project_fastqs/asap_fastqs/$sample" -c "$threads"

                mkdir -p temp

                kallisto bus -i adt_index/FeaturesMismatch.idx -o temp/ -x 0,0,16:0,16,26:1,0,0 -t "$threads" "$project_fastqs/asap_fastqs/$sample"*

                bustools correct -w "$atac_whitelist" temp/output.bus -o temp/output_corrected.bus
                bustools sort -t "$threads" -o temp/output_sorted.bus temp/output_corrected.bus
                bustools count -o "$sample/adt/" --genecounts -g adt_index/FeaturesMismatch.t2g -e temp/matrix.ec -t temp/transcripts.txt temp/output_sorted.bus

                rm -r temp
        } 2>&1
done
```
