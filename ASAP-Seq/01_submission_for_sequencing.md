# Submitting prepared libraries for sequencing

[Two types of libraries](https://github.com/romagnanilab/single_cell_analysis/blob/2fc8ac359c6e9da703de873a50b26d6650b4b64d/static_files/asap_original_protocol.pdf) are generated following the [ASAP-seq](https://cite-seq.com/asapseq/) experimental pipeline:
1. One for the peaks, with the library structure [here](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html)
2. Another for the antibody capture and hashtag information which you can find the library strcuture of [here](https://github.com/romagnanilab/single_cell_analysis/blob/3c8e77babbde4d7305c3c9f8d893ef4104df45ed/static_files/asap_assay_scheme.pdf).

ATAC-seq libaries generally require less sequencing saturation that RNA-seq libraries due to their decreased complexity, deriving from their higher signal-to-noise ratio and sparse coverage across the large genome. 10X recommend to sequence to [25,000 read pairs per nucleus](https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/sequencing/sequencing-requirements-for-single-cell-atac), but we have generally gone higher at 50,000 or so. 

When you submit your libraries to get sequenced, take care in choosing the read lengths for Read 1, the i7 and i5 indices, and Read 2.
