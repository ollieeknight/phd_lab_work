# Submitting prepared libraries for sequencing


[Two types of libraries](https://github.com/romagnanilab/single_cell_analysis/blob/2fc8ac359c6e9da703de873a50b26d6650b4b64d/static_files/asap_original_protocol.pdf) are generated following the [ASAP-seq](https://cite-seq.com/asapseq/) experimental pipeline:
1. One for the peaks, with the library structure [here](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_scATAC.html)
2. Another for the antibody capture and hashtag information which you can find the library strcuture of [here](https://github.com/romagnanilab/single_cell_analysis/blob/3c8e77babbde4d7305c3c9f8d893ef4104df45ed/static_files/asap_assay_scheme.pdf).

ATAC-seq libaries generally require less sequencing saturation that RNA-seq libraries due to their decreased complexity, deriving from their higher signal-to-noise ratio and sparse coverage across the large genome. 10X recommend to sequence to [25,000 read pairs per nucleus](https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/sequencing/sequencing-requirements-for-single-cell-atac), but we have generally gone higher to 50,000 or so.

For the antibody capture, the quality of data is generally lower due to the alternative method of how we capture antibody reads; using the bridge oligo versus the read being captured inside the GEM like it is with scRNA-seq and multiome protocols. ** ADD MORE HERE **

When you submit your libraries to get sequenced, take care in choosing the read lengths for Read 1, the i7 and i5 indices, and Read 2. You can find these [here](https://www.10xgenomics.com/support/single-cell-atac/documentation/steps/sequencing/sequencing-requirements-for-single-cell-atac), but the bottom line is to sequence *at least*:

| Read | Length |
| ------------- | ------------- |
| Read 1 | 50 |
| i7 index | 8 |
| i5 index | 16  |
| Read 2 | 50  |

Using an SP or S1 flowcell, we often sequence longer, choosing 100, 8, 16, 50. This will be important later while using `cellranger`'s `bcl2fastq` wrapper and  `mgatk` for mitochondrial genotyping. 

## The process of sequencing

Here's a general overview of sequencing, which might help you conceptualise how your libraries become files for us to process in the next step:

1. Loading a flowcell: A flowcell is a glass slide or chip that contains millions of small wells or channels where DNA fragments will be sequenced. The library is loaded onto the flowcell and the adapters of each cDNA fragment will bind to the flowcell surface.

2. Cluster generation: fragments are amplified and spatially localised on the flowcell surface. This amplification process involves the addition of enzymes, primers, and nucleotides which allow DNA replication. Clusters of identical DNA fragments are formed in each well or channel of the flowcell.

3. Sequencing-by-Synthesis (SBS): During illumina sequencing, fluorescently labeled nucleotides are added to the flowcell, and the sequencing instrument detects the fluorescence emitted by each incorporated nucleotide. The nucleotides are added one at a time, and their fluorescence signals are recorded. This process allows the determination of the DNA sequence.

4. Image acquisition: During the SBS process, a series of images are captured by the sequencing instrument. These images correspond to the fluorescent signals emitted by the nucleotides as they are incorporated into the growing DNA strands. The images capture the identity and location of each nucleotide in the clusters.

5. Base calling: Base calling is the process of translating the fluorescent signals from the images into nucleotide sequence information. Software analyses the images and assigns a base call (A, C, G, T) to each location in the clusters based on the observed fluorescence patterns. Quality scores are also assigned to each base call, reflecting the confidence in the accuracy of the base call at that specific position.

6. BCL file generation: Once the base calling and quality scoring are completed for all clusters, the information is compiled into a BCL file.
