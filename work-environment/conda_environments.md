# Setting up conda environments for processing genomics files
Once you've logged into the BIH-CUBI cluster through the command line, setting up your conda environments is crucial for processing genetic data.
```shell
cd /fast/work/users/${USER}/ && mkdir bin/ && cd bin/

# download, install, and update miniconda 
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p miniconda3 && rm Miniconda3-latest-Linux-x86_64.sh
source miniconda3/etc/profile.d/conda.sh && conda activate

# modify where conda looks for packages
nano ~/.condarc

# copy and paste this into nano text editor
channels:
  - conda-forge
  - bioconda
  - defaults
show_channel_urls: true
changeps1: true
channel_priority: strict
# CTRL + X, then Y, then enter to leave

conda upgrade --all 
y

nano ~/.bashrc
#paste this: source /fast/work/users/$USER/bin/miniconda3/etc/profile.d/conda.sh
# CTRL + X, then Y, then enter to leave

```

## Processing raw BCL files
```shell
conda create -y -n bcl_to_fastq -c bih-cubi bcl2fastq2
```

## Genotyping genomic BAM files
```shell
conda create -y -n donor_genotyping python cellsnp-lite  
conda activate donor_genotyping
pip install vireoSNP
```

## Doublet detection from scATAC fragment overlap
```shell
conda create -y -n amulet_overlap numpy=1.19 pandas scipy statsmodels
```

## Genotyping scATAC mitochondrial DNA
```shell
conda create -y -n mitochondrial_genotyping openjdk r-base=4.2.3 r-data.table r-matrix bioconductor-genomicranges bioconductor-summarizedexperiment
conda activate mitochondrial_genotyping
pip install mgatk
```

## Collating antibody capture counts per cell
```shell
conda create -y -n adt_count python kallisto bustools 
conda activate adt_count
pip install bio
```

## Cellbender
```shell
conda create -y -n cellbender -c nvidia python=3.7 cuda-toolkit cuda-nvcc
conda activate cellbender
pip install cellbender
```

## Manipulating genomic files

```shell
conda create -y -n genome_processing bcftools samtools bedtools bwa
```

## Using python packages in `R` through `reticulate`

```shell
conda create -y -n R_sc_python numpy leidenalg umap-learn macs2 scanpy scvi-tools
```

Then, in R, start your script with
```R
Sys.setenv(RETICULATE_MINICONDA_PATH = '~/bin/miniconda3/')
Sys.setenv(PATH = paste('~/bin/miniconda3/envs/R_sc_python/lib/python3.11/site-packages/', Sys.getenv()['PATH'], sep = ':'))
library(reticulate)
use_miniconda('~/bin/miniconda3/envs/R_sc_python')
```

And for MACS2 peaks calling:
```R
peaks <- Signac::CallPeaks(alldata, assay = 'ATAC', macs2.path = '~/bin/miniconda3/envs/R_sc_python/bin/macs2')
```
Where `alldata` is your seurat object.
