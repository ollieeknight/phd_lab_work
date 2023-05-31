# Setting up your working directory
How you choose to organise your files is entirely your choice, but I've provided here a starting guide which I use for mine.
```
$ tree ~/
├── group -> /fast/home/groups/ag_romagnani/
    ├── scratch -> /fast/scratch/groups/ag_romagnani
    └── work -> /fast/work/groups/ag_romagnani
├── ondemand -> work/bin/ondemand
├── scratch -> /fast/scratch/users/knighto_c
    ├── BIH_TRASH
    ├── ngs
    ├── nf_work
    └── tmp
└── work -> /fast/work/users/knighto_c/
    ├── bin
    ├── data
    └── slurm

5 directories, 0 files
```
You will of course notice that the immediate folders are symbolic links (created with `ln -s actual_folder target_folder`) to other storage locations present on the cluster:
1. You start in `home`, which is 1Gb in size
2. The group `home` folder was created with `ln -s /fast/home/groups/ag_romagnani/ group`, for quick access
3. The `ondemand` folder was linked with `ln -s ~/work/bin/ondemand ondemand` group, for the `ondemand` program to find
4. The `scratch` folder, for processing large genetic datasets, with a 200Tb file and 2.2M file number size limit
5. The `work` folder, for smaller projects, with a 1Tb file size limit

## scratch
### BIH_TRASH
This folder will automatically be present and is run by the CUBI admin - any file with an age older than 2 weeks is moved here, and be permanently removed in another two weeks. This is to maintain open space on the `scratch` storage space.

### nf-work
You can create this folder as per the tutorial [here](https://github.com/ollieeknight/single_cell_analysis/blob/main/work-environment/nextflow.md), for the `nextflow` program.

### tmp
This is an important folder which programs will use for temporary write storage, which you should add to your login `.bashrc` file with:
```
export TMPDIR=/fast/users/$USER/scratch/tmp
mkdir -p $TMPDIR/nf_work
```

### ngs
You can create this folder yourself with `mkdir -p ~/scratch/ngs`, as an isolated space to download and process genetic data. 

## work
### bin
Within my bin folder I have [`miniconda3`](https://github.com/ollieeknight/single_cell_analysis/blob/main/work-environment/conda_environments.md), `nextflow`, and `ondemand` installed. These are folders which are for my personal programs. Shared static programs like `cellranger` are found in `/fast/home/groups/ag_romagnani/bin/`.

### data
Within my data folder I have my project files, for example `NK_CMVpos` with `data`, `objects`, and `scripts` folders. The data folder has **no** `bam` files in it as these take up a lot of space.

### slurm
This is a folder where I set `slurm` log and error messages to output.
