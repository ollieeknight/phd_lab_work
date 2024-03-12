#!/bin/bash

echo ""
echo "Setting up bedtools env"
echo ""

# Conda environment setup
conda create -y -p "$TMPDIR/reference" bedtools || { echo "Error creating conda environment; is miniconda in your path?"; exit 1; }

echo ""

# Check cellranger-arc path
cellranger_arc_path="$HOME/group/work/bin/cellranger-arc-2.0.2"
[ ! -d "$cellranger_arc_path" ] || { echo "Error: cellranger-arc path not found"; exit 1; }
export PATH="$cellranger_arc_path:$PATH"

echo "Do you want to create the reference genome here ($PWD)? (y/newdir/n)"
read -r choice

# Validate user input
while [[ ! $choice =~ ^[YyNn]$ && $choice != "newdir" ]]; do
    echo "Invalid input. Please enter Y, N, or newdir"
    read -r choice
done

if [[ $choice == "newdir" ]]; then
    echo ""
    echo "Enter the directory:"
    read -r ref_directory

    if [ ! -d "$ref_directory" ]; then
        echo "Error: Directory not found, please create it"
        exit 1
    fi
fi

if [[ $choice =~ ^[Yy]$ ]]; then
    echo ""
    echo "Downloading reference files, this might take a while..."
    echo ""
    # Add commands for downloading reference files
elif [[ $choice =~ ^[Nn]$ ]]; then
    echo ""
    echo "Not downloading reference files, exiting"
    echo ""
    exit
fi

cd "$ref_directory" || { echo "Error: Could not change to the specified directory"; exit 1; }

# Genome setup
genome="GRCh38-hardmasked-optimised-arc"
build="${genome}-build"
source="${genome}-source"
mkdir -p "${build}" "${source}"

# Download files
download_file() {
    local url=$1
    local file=$2
    [ -f "$file" ] || curl -sS "$url" > "$file"
}

fasta_name="Homo_sapiens.GRCh38.dna_sm.primary_assembly"
download_file "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/${fasta_name}.fa.gz" "${source}/${fasta_name}.fa"
download_file "https://storage.googleapis.com/generecovery/human_GRCh38_optimized_annotation_v2.gtf.gz" "${source}/human_GRCh38_optimized_annotation_v2.gtf"
download_file "https://testjaspar.uio.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt" "${source}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt"

# Modify fasta and motifs
fasta_mod="${build}/${fasta_name}_hardmasked.fa.mod"
cat "${source}/${fasta_name}.fa" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "${fasta_mod}"

motifs_mod="${build}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt.mod"
awk '{print ">" $2 "_" substr($1,2)}' "${source}/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt" > "${motifs_mod}"

# Mask fasta
conda run -p $TMPDIR/reference bedtools maskfasta -fi "${fasta_mod}" -bed "${source}/hg38.full.blacklist.bed" -fo "${fasta_mod}"

# Genome configuration
config_in="${build}/genome.config"
cat > "${config_in}" <<EOF
{
    organism: "Homo_sapiens",
    genome: ["${genome}"],
    input_fasta: ["${fasta_mod}"],
    input_gtf: ["${source}/human_GRCh38_optimized_annotation_v2.gtf"],
    input_motifs: "${motifs_mod}",
    non_nuclear_contigs: ["chrM"]
}
EOF

# Cellranger-arc mkref
cellranger-arc mkref --ref-version 'A' --config "${config_in}"

# Cleanup
rm -r "${source}" "${build}"

