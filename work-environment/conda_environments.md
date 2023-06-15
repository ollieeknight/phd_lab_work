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

conda install mamba
```

## Processing raw BCL files
```shell
mamba create -y -n bcl_to_fastq -c bih-cubi bcl2fastq2
```

## Genotyping genomic BAM files
```shell
mamba create -y -n donor_genotyping python cellsnp-lite  
conda activate donor_genotyping
pip install vireoSNP
```

## Doublet detection from scATAC fragment overlap
```shell
mamba create -y -n amulet_overlap numpy pandas scipy statsmodels
```

## Genotyping scATAC mitochondrial DNA
```shell
mamba create -y -n mitochondrial_genotyping openjdk r-data.table r-matrix bioconductor-genomicranges bioconductor-summarizedexperiment
conda activate mitochondrial_genotyping
pip install mgatk
```
## Collating antibody capture counts per cell
```shell
mamba create -y -n adt_count python kallisto bustools 
conda activate adt_count
pip install bio
```

## Cellbender
```shell
mamba create -y -n cellbender -c nvidia python=3.7 cuda-toolkit cuda-nvcc
conda activate cellbender
pip install git+https://github.com/ollieeknight/CellBender@sf_dev_0.3.0_postreg
```

## Manipulating vcf files

```shell
mamba create -y -n vcf_process bcftools 
```

## Using python packages in `R` through `reticulate`

```shell
mamba create -y -n r_reticulate_python numpy leidenalg umap-learn macs2 scanpy scvi-tools
conda activate sc_python
```

Then, in R, start your script with
```R
Sys.setenv(RETICULATE_MINICONDA_PATH = '/fast/work/users/$username/bin/miniconda3/')
Sys.setenv(PATH = paste('/fast/work/users/$username/bin/miniconda3/envs/r_reticulate_python/lib/python3.11/site-packages/', Sys.getenv()['PATH'], sep = ':'))
library(reticulate)
use_miniconda('/fast/work/users/$username/bin/miniconda3/envs/r_reticulate_python')
```
Replacing `$username` with your username.

And for MACS2 peaks calling:
```R
peaks <- Signac::CallPeaks(alldata, assay = 'ATAC', macs2.path = '/fast/work/users/$username/bin/miniconda3/envs/peaks/bin/macs2')
```
Where `alldata` is your seurat object, and replacing `$username`, again.
