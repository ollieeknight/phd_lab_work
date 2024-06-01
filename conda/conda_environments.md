# Setting up conda environments for processing genomics files
Once you've logged into the BIH-CUBI cluster through the command line, setting up your conda environments is crucial for processing genetic data.
```shell
mkdir $HOME/work/bin/ && cd $HOME/work/bin/

# download, install, and update miniconda 
curl -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/work/bin/miniconda3 && rm Miniconda3-latest-Linux-x86_64.sh
source miniconda3/etc/profile.d/conda.sh && conda activate

# modify conda repositories  
nano ~/.condarc

# copy and paste this into nano (CTRL+C here, right click to paste)
channels:
  - conda-forge
  - bioconda
  - defaults
show_channel_urls: true
changeps1: true
channel_priority: flexible
# close by CTRL+X and y and enter

conda upgrade --all -y
conda config --set solver libmamba

```

## Processing raw BCL files
```bash
conda create -y -n bcl2fastq2 -c bih-cubi bcl2fastq2
```

## Genotyping genomic BAM files
```bash
conda create -y -n vireo python cellsnp-lite  
conda activate vireo
pip install vireoSNP
```

## Doublet detection from scATAC fragment overlap
```bash
conda create -y -n amulet numpy=1.19 pandas scipy statsmodels

# Then install AMULET via 
git clone https://github.com/UcarLab/AMULET
```

## Genotyping scATAC mitochondrial DNA
```bash
conda create -y -n mgatk openjdk r-base=4.2.3 r-data.table r-matrix bioconductor-genomicranges bioconductor-summarizedexperiment
conda activate mgatk
pip install mgatk

# Important to note: mgatk requires ~8gb RAM per core and so it's best to run with -c 8 and >64 GB RAM 
```

## Collating antibody capture counts per cell
```bash
conda create -y -n kite python kallisto bustools 
conda activate kite
pip install bio
```

## Cellbender
```bash
conda create -y -n cellbender -c nvidia python=3.7 cuda-toolkit cuda-nvcc
conda activate cellbender
pip install cellbender
```

## Manipulating genomic files

```bash
conda create -y -n vcf bcftools samtools bedtools bwa
```

## Using python packages in `R` through `reticulate`

```bash
conda create -y -n r-reticulate -c vtraag python-igraph pandas umap-learn scanpy macs2 scvi-tools
conda activate r-reticulate
conda install -c vtraag leidenalg
```

Then, in R, start your script with
```R
Sys.setenv(RETICULATE_MINICONDA_PATH = '~/work/bin/miniconda3/')
Sys.setenv(PATH = paste('~/work/bin/miniconda3/envs/r-reticulate/lib/python3.10/site-packages/', Sys.getenv()['PATH'], sep = ':'))
library(reticulate)
use_miniconda('~/work/bin/miniconda3/envs/r-reticulate/')
```

And for MACS2 peaks calling:
```R
peaks <- Signac::CallPeaks(alldata, assay = 'ATAC', macs2.path = '~/bin/miniconda3/envs/r-reticulate/bin/macs2')
```
Where `alldata` is your seurat object, with your 'ATAC' (peaks) object set as the default assay.
