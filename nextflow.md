# Using nextflow for processing WES data

## Install nextflow and dependencies

```bash
# Install Java
mkdir -p "$HOME/work/bin"
cd "$HOME/work/bin" || exit
wget https://download.oracle.com/java/17/archive/jdk-17.0.7_linux-x64_bin.tar.gz
tar -xf jdk-17.0.7_linux-x64_bin.tar.gz && rm jdk-17.0.7_linux-x64_bin.tar.gz
export PATH="$HOME/work/bin/jdk-17.0.7/bin:$PATH"

# Install Nextflow
mkdir -p "$HOME/work/bin/nextflow"
cd "$HOME/work/bin/nextflow" || exit
wget -qO- https://get.nextflow.io | bash
export PATH="$HOME/work/bin/nextflow:$PATH"
nextflow self-update

# Setup directories
mkdir -p "$HOME/scratch/tmp/nextflow"
mkdir -p "$HOME/work/bin/nextflow/apptainer_cache"
```

```bash
cat <<EOL >> ~/.bashrc
# Nextflow environment variables
export NXF_HOME=${HOME}/work/bin/nextflow
export NXF_JAVA_HOME=${HOME}/work/bin/jdk-17.0.7
export NXF_TEMP=${HOME}/scratch/tmp/nextflow
export NXF_APPTAINER_CACHEDIR=${HOME}/work/bin/nextflow/apptainer_cache
export NXF_WORK=${HOME}/scratch/tmp/nextflow
export APPTAINERENV_NXF_TASK_WORKDIR=${HOME}/scratch/tmp/nextflow
export APPTAINERENV_TMPDIR=${HOME}/scratch/tmp/nextflow
export NXF_APPTAINER_HOME_MOUNT=true
EOL
```

### Configuration file - you can also find this under ~/group/ref/WES/

```bash
executor {
  name = 'slurm'
  queueSize = 32
}

profiles {
    apptainer {
        conda.enabled           = false
        apptainer.enabled       = true
        apptainer.autoMounts    = true
        docker.enabled          = false
        apptainer.runOptions    = '-B /fast'
    }
}

env {
    TMPDIR = '/fast/scratch/users/knighto_c/tmp/nextflow'
    APPTAINER_CACHEDIR = '/fast/work/users/knighto_c/bin/nextflow/apptainer_cache'
    APPTAINERENV_NXF_TASK_WORKDIR = '/fast/scratch/users/knighto_c/tmp/nextflow'
    APPTAINERENV_TMPDIR = '/fast/scratch/users/knighto_c/tmp/nextflow'
}

params {
    max_memory = 16.GB
    max_cpus = 8
    max_time = 48.h
    genomes {
       'GATK.GRCh38' {
            ascat_alleles         = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_alleles_hg38.zip"
            ascat_genome          = 'hg38'
            ascat_loci            = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_loci_hg38.zip"
            ascat_loci_gc         = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/GC_G1000_hg38.zip"
            ascat_loci_rt         = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/RT_G1000_hg38.zip"
            bwa                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/"
            bwamem2               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index/"
            dragmap               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/dragmap/"
            chr_dir               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes"
            cf_chrom_len          = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/Length/Homo_sapiens_assembly38.len"
            dbsnp                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"
            dbsnp_tbi             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz.tbi"
            dbsnp_vqsr            = '--resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz'
            dict                  = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict"
            fasta                 = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
            fasta_fai             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai"
            germline_resource     = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz"
            germline_resource_tbi = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz.tbi"
            intervals             = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed"
            known_snps            = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz"
            known_snps_tbi        = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi"
            known_snps_vqsr       = '--resource:1000G,known=false,training=true,truth=true,prior=10.0 1000G_omni2.5.hg38.vcf.gz'
            known_indels          = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz"
            known_indels_tbi      = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz.tbi"
            known_indels_vqsr     = '--resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
            mappability           = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/Control-FREEC/out100m2_hg38.gem"
            pon                   = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz"
            pon_tbi               = "${params.igenomes_base}/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi"
            snpeff_db             = 105
            snpeff_genome         = 'GRCh38'
            snpeff_version        = '5.1'
            vep_cache_version     = 106
            vep_genome            = 'GRCh38'
            vep_species           = 'homo_sapiens'
            vep_version           = '106.1'
        }
   }
}
```

### Sample sheet

```csv
patient,status,sample,lane,fastq_1,fastq_2
PATIENT1,0,CTR01,1,CTR01_R1.fastq.gz,CTR01_R2.fastq.gz
PATIENT1,1,TMR01,1,TMR01_R1.fastq.gz.TMR01_R1.fastq.gz
```

### Run command
```
nextflow run nf-core/sarek -r 3.4.0 --input samplesheets/hc_paired_samplesheet.csv --outdir outs/HC_paired_outs --genome GATK.GRCh38 --igenomes_base ~/group/work/ref/igenomes -profile apptainer -c ~/group/work/ref/WES/cubi.config --tools mutect2,strelka,vep --intervals ~/group/work/ref/WES/hg38_exome_v2.0.2_targets_sorted_validated.re_annotated.bed --wes --only_paired_variant_calling
```

## Extracting tables from vcf files
One you have your annotated variants from the `sarek` pipeline, we can extract information from `mutect2`- and `strelka`-run paired samples using bcftools

```bash
conda create -y -n vcf bcftools samtools bedtools bwa
conda activate
```

### Strelka

```bash
sample=xyz
echo -e "CHROM\tPOS\tREF\tALT\tA\tC\tG\tT\tFILTER\t$(bcftools +split-vep -l *snvs_VEP.ann.vcf.gz | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > paired_${sample}.tsv && bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t[_%AU]\t[_%CU]\t[_%GU]\t[_%TU]\t%FILTER\t%CSQ\n' -d -A tab *snvs_VEP.ann.vcf.gz >> paired_${sample_}.tsv
```

This is agnostic of how you label your paired samples - \*snvs_VEP.ann.vcf.gz wildcards the names. 

### Mutect2

```bash
sample=xyz
echo -e "CHROM\tPOS\tREF\tALT\tAF\tDP\tFILTER\t$(bcftools +split-vep -l *_VEP.ann.vcf.gz | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > paired_${sample}.tsv && bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t[_%AF]\t[_%DP]\t%FILTER\t%CSQ\n' -d -A tab *_VEP.ann.vcf.gz >> paired_${sample}.tsv
```

This is agnostic of how you label your paired samples - \*_VEP.ann.vcf.gz wildcards the names. 
