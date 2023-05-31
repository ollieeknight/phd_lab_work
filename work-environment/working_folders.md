# Setting up your working directory
How you choose to organise your files is entirely your choice, but I've provided here a starting guide which I use for mine.
```
$ tree ~/
├── group -> /fast/home/groups/ag_romagnani/
    ├── scratch -> /fast/scratch/groups/ag_romagnani
    └── work -> /fast/work/groups/ag_romagnani
├── ondemand -> work/bin/ondemand
├── R -> work/bin/miniconda3/envs/sc_R/lib/R/library/
├── scratch -> /fast/scratch/users/knighto_c
    ├── BIH_TRASH
    ├── ngs
    └── tmp
└── work -> /fast/work/users/knighto_c/
    ├── bin
    ├── data
    └── slurm

5 directories, 0 files
```
Directories which have arrows pointing to another are symbolic links, created with the command `ln -s actual_folder link_folder`

## scratch

## work
### bin
Within my bin folder I have `miniconda3`, `nextflow`, and `ondemand` installed. These are folders which are for my personal programs. Shared static programs like `cellranger` are found in `/fast/home/groups/ag_romagnani/bin/`.

### data
Within my data folder I have my project files, for example `NK_CMVpos` with `data`, `objects`, and `scripts` folders. The data folder has **no** `bam` files in it as these take up a lot of space.

### slurm
This is a folder where I set `slurm` log and error messages to output.
