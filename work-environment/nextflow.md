# Using nextflow for high-throughput genetic data processing
For now I am only going to put my environment variables here. Talk to me if you would like any more guidance!

## Install nextflow and modify dependencies

```
cd ~/work/bin
wget https://download.oracle.com/java/17/archive/jdk-17.0.7_linux-x64_bin.tar.gz 
tar -xf jdk-17.0.7_linux-x64_bin.tar.gz && rm jdk-17.0.7_linux-x64_bin.tar.gz
export PATH=/data/gpfs-1/users/knighto_c/work/bin/jdk-17.0.7/bin:$PATH
mkdir -p nextflow
cd nextflow
wget -qO- https://get.nextflow.io | bash
export PATH=/fast/home/users/knighto_c/work/bin/nextflow:$PATH

mkdir -p /fast/home/users/knighto_c/scratch/nf_work
mkdir -p /fast/home/users/knighto_c/work/bin/nextflow/singularity_cache

export NXF_HOME=/data/gpfs-1/users/knighto_c/work/bin/nextflow
export NXF_SINGULARITY_CACHEDIR=/data/gpfs-1/users/knighto_c/work/bin/nextflow/singularity_cache
export NXF_WORK=/fast/users/knighto_c/scratch/tmp/nf_work
export NXF_JAVA_HOME=/data/gpfs-1/users/knighto_c/work/bin/jdk-17.0.7
```

## Processing Whole Exome Sequencing (WES) data with [Sarek](https://nf-co.re/sarek)

### Configuration file

```
executor {
  name = 'slurm'
  queueSize = 16
}
singularity {
  enabled = true
}
env {
    TMPDIR = '/fast/scratch/users/knighto_c/tmp/nf_work'
    SINGULARITY_CACHEDIR = '/fast/work/users/knighto_c/bin/nextflow/singularity_cache'
}
params {
    max_memory = 64.GB
    max_cpus = 32
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
            snpeff_db             = 'GRCh38.105'
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
```
patient,status,sample,lane,fastq_1,fastq_2
LGL062,0,CT062,1,/data/gpfs-1/users/knighto_c/work/data/S8163/S8163_fastq/S000021_S8163Nr6.1.fastq.gz,/data/gpfs-1/users/knighto_c/work/data/S8163/S8163_fastq/S000021_S8163Nr6.2.fastq.gz
LGL062,1,NK062,1,/data/gpfs-1/users/knighto_c/work/data/S8163/S8163_fastq/S000021_S8163Nr5.1.fastq.gz,/data/gpfs-1/users/knighto_c/work/data/S8163/S8163_fastq/S000021_S8163Nr5.2.fastq.gz
```
### Run command
```
nextflow run nf-core/sarek -r 3.1.2 --input samplesheet.csv --outdir project_id_outs --genome GATK.GRCh38 -profile singularity -c cubi.config --igenomes_base /fast/home/groups/ag_romagnani/work/ref/igenomes --tools mutect2,vep
```

## Extracting tables from vcf files
One you have your annotated variants from the `sarek` pipeline, we can extract information from `mutect2`- and `strelka`-run paired samples. You will need to activate your vcf conda environment, which we created [here](https://github.com/ollieeknight/single_cell_analysis/blob/main/work-environment/conda_environments.md#manipulating-vcf-files) 

### Strelka
```
sample=xyz
echo -e "CHROM\tPOS\tREF\tALT\tA\tC\tG\tT\tFILTER\t$(bcftools +split-vep -l *snvs_VEP.ann.vcf.gz | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > paired_${sample}.tsv && bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t[_%AU]\t[_%CU]\t[_%GU]\t[_%TU]\t%FILTER\t%CSQ\n' -d -A tab *snvs_VEP.ann.vcf.gz >> paired_${sample_}.tsv
```
This is agnostic of how you label your paired samples - \*snvs_VEP.ann.vcf.gz wildcards the names. 

### Mutect2
```
sample=xyz
echo -e "CHROM\tPOS\tREF\tALT\tAF\tDP\tFILTER\t$(bcftools +split-vep -l *_VEP.ann.vcf.gz | cut -f 2 | tr '\n' '\t' | sed 's/\t$//')" > paired_${sample}.tsv && bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t[_%AF]\t[_%DP]\t%FILTER\t%CSQ\n' -d -A tab *_VEP.ann.vcf.gz >> paired_${sample}.tsv
```
This is agnostic of how you label your paired samples - \*_VEP.ann.vcf.gz wildcards the names. 


