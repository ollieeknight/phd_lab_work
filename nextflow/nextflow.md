# Running Nextflow Sarek for WGS/WES data

- [nf-co.re/sarek/latest/](https://nf-co.re/sarek/latest/)
- [github.com/nf-core/sarek](https://github.com/nf-core/sarek)

## Install Nextflow and dependencies

### Setup directories

```bash
mkdir -p "$HOME/work/bin"
mkdir -p "$HOME/work/bin/nextflow"
mkdir -p "$HOME/scratch/tmp/nextflow"
mkdir -p "$HOME/work/bin/nextflow/apptainer_cache"
```

### Install Java

```bash
cd "$HOME/work/bin" || exit
wget https://download.oracle.com/java/21/archive/jdk-21.0.3_linux-x64_bin.tar.gz
tar -xf jdk-21.0.3_linux-x64_bin.tar.gz && rm jdk-21.0.3_linux-x64_bin.tar.gz
export PATH="$HOME/work/bin/jdk-21.0.3/bin:$PATH"
```

### Install Nextflow

```bash
cd "$HOME/work/bin/nextflow" || exit
wget -qO- https://get.nextflow.io | bash
export PATH="$HOME/work/bin/nextflow:$PATH"
nextflow self-update
```


## Set up Nextflow environment variables

```bash
cat <<EOL >> ~/.bashrc
export NXF_HOME=${HOME}/work/bin/nextflow
export NXF_JAVA_HOME=${HOME}/work/bin/jdk-21.0.3
export NXF_TEMP=${HOME}/scratch/tmp/nextflow
export NXF_APPTAINER_CACHEDIR=${HOME}/work/bin/nextflow/apptainer_cache
export NXF_WORK=${HOME}/scratch/tmp/nextflow
export APPTAINERENV_NXF_TASK_WORKDIR=${HOME}/scratch/tmp/nextflow
export APPTAINERENV_TMPDIR=${HOME}/scratch/tmp/nextflow
export NXF_APPTAINER_HOME_MOUNT=true
EOL

source ~/.bashrc
```

### Configuration file - take from `cubi_wes.config` or `cubi_wgs.config`, depending on your use case!

- [Whole exome sequencing config file](https://github.com/ollieeknight/phd_lab_work/blob/main/nextflow/cubi_wes.config)
- [Whole genome sequencing config file](https://github.com/ollieeknight/phd_lab_work/blob/main/nextflow/cubi_wgs.config)

## Running Nextflow sarek

### Example samplesheet.csv

```csv
patient,status,sample,lane,fastq_1,fastq_2
PATIENT1,0,CTR01,1,CTR01_R1.fastq.gz,CTR01_R2.fastq.gz
PATIENT1,1,TMR01,1,TMR01_R1.fastq.gz.TMR01_R1.fastq.gz
```

### Whole exome sequencing

#### Germline

```bash
NXF_VER=24.04.4 nextflow run nf-core/sarek -r 3.4.4 \
  --input samplesheet.csv \
  --outdir outs \
  --genome GATK.GRCh38 \
  -profile apptainer \
  -c /data/cephfs-1/home/users/knighto_c/group/work/ref/wes/cubi_wes.config \
  --tools deepvariant,cnvkit,vep \
  --wes
```

#### Tumour

```bash
NXF_VER=24.04.4 nextflow run nf-core/sarek -r 3.4.4 \
  --input samplesheet.csv \
  --outdir outs \
  --genome GATK.GRCh38 \
  -profile apptainer \
  -c /data/cephfs-1/home/users/knighto_c/group/work/ref/wes/cubi_wes.config \
  --tools strelka,mutect2,vep \
  --wes \
  --only_paired_variant_calling
```

### Whole genome sequencing

```bash
NXF_VER=24.04.4 nextflow run nf-core/sarek -r 3.4.4 \
  --input samplesheet.csv \
  --outdir outs \
  --genome GATK.GRCh38 \
  -profile apptainer \
  -c /data/cephfs-1/home/users/knighto_c/group/work/ref/wes/cubi_wgs.config \
  --tools strelka,mutect2,vep \
  --only_paired_variant_calling
```
