**s**```bash
conda create -y -n genome_processing bcftools samtools bedtools bwa
conda activate genome_processing

export PATH=~/group/work/bin/cellranger-arc-2.0.2:$PATH

genome="GRCh38-hardmasked-optimised" # enter end folder name here
build="${genome}-build"
mkdir -p "$build"

source=${genome}-source
mkdir -p $source
fasta_url="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz"
fasta_in="$source/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf_url="https://storage.googleapis.com/generecovery/human_GRCh38_optimized_annotation_v2.gtf.gz"
gtf_in="${source}/human_GRCh38_optimized_annotation_v2.gtf"
motifs_url="https://testjaspar.uio.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"
motifs_in="${source}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi

if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi

if [ ! -f "$motifs_in" ]; then
    curl -sS "$motifs_url" > "$motifs_in"
fi

fasta_mod="$build/$(basename "$fasta_in").mod"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_mod"

motifs_mod="$build/$(basename "$motifs_in").mod"
awk '{
    if ( substr($1, 1, 1) == ">" ) {
        print ">" $2 "_" substr($1,2)
    } else {
        print
    }
}' "$motifs_in" > "$motifs_mod"

wget https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed -o $source/hg38.full.blacklist.bed
mv $build/Homo_sapiens.GRCh38.dna.primary_assembly.fa.mod $build/Homo_sapiens.GRCh38.dna.primary_assembly_original.fa.mod
bedtools maskfasta -fi $build/Homo_sapiens.GRCh38.dna.primary_assembly_original.fa.mod -bed $source/hg38.full.blacklist.bed -fo $build/Homo_sapiens.GRCh38.dna.primary_assembly_hardmasked.fa.mod

config_in="$build/genome.config"
echo """{
    organism: \"Homo_sapiens\"
    genome: [\""$genome"\"]
    input_fasta: [\""$build/Homo_sapiens.GRCh38.dna.primary_assembly_hardmasked.fa.mod"\"]
    input_gtf: [\""$gtf_in"\"]
    input_motifs: \""$motifs_mod"\"
    non_nuclear_contigs: [\"chrM\"]
}""" > "$config_in"

cellranger-arc mkref --ref-version 'A' --config $config_in

rm -r $source $build
```
