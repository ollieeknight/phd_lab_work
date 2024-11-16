pyscenic grn RA_Immune_filtered.loom --num_workers 64 -o adj.tsv ~/work/bin/scenic/allTFs_hg38.txt

f_db_glob="/data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/scenic/*feather"
f_db_names=$(echo $f_db_glob)

pyscenic ctx adj.tsv $f_db_names --annotations_fname /data/cephfs-1/work/groups/romagnani/users/knighto_c/bin/scenic/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname RA_Immune_filtered.loom --output reg.csv --no_pruning --mask_dropouts --num_workers 64

pyscenic aucell RA_Immune_filtered.loom reg.csv --output RA_Immune_filtered_scenic.loom --num_workers 64 --auc_threshold=0.01