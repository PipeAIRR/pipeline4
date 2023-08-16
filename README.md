The pipeline process BCR repertoire sequence that were produce in the same fashion as those in [Galson et al. 2020](https://www.frontiersin.org/articles/10.3389/fimmu.2020.605170)

Library preperation and sequencing method:


RNA Prep and BCR Sequencing

Total RNA from 5 × 106 PBMCs was isolated using RNeasy kits (Qiagen). First-strand cDNA was generated from total RNA using SuperScript RT IV (Invitrogen) and IgA and IgG isotype specific primers (17) including UMIs at 50°C for 45 min (inactivation at 80°C for 10 min).
The resulting cDNA was used as template for High Fidelity PCR amplification (KAPA, Roche) using a set of 6 FR1-specific forward primers (17) including sample-specific barcode sequences (6 bp) and a reverse primer specific to the RT primer (initial denaturation at 95°C for 3 min, 25 cycles at 98°C for 20 s, 60°C for 30 s, 72°C for 1 min and final extension at 72°C for 7 min). The amount of BCR heavy chain amplicons (~450 bp) was quantified by TapeStation (Beckman Coulter) and gel-purified.
Dual-indexed sequencing adapters (KAPA) were ligated onto 500-ng amplicons per patient using the HyperPrep library construction kit (KAPA) and the adapter-ligated libraries were finally PCR-amplified for 3 cycles (98°C for 15 s, 60°C for 30 s, 72°C for 30s, final extension at 72°C for 1 min). Pools of 10, 9, and 12 libraries were sequenced across three runs on an Illumina MiSeq using 2 × 300 bp chemistry.

The Immcantation framework (docker container v3.0.0) was used for sequence processing (18, 19). Briefly, paired-end reads were joined based on a minimum overlap of 20 nt, and a max error of 0.2, and reads with a mean phred score below 20 were removed. Primer regions, including UMIs and sample barcodes, were then identified within each read, and trimmed. Together, the sample barcode, UMI, and constant region primer were used to assign molecular groupings for each read. Within each grouping, usearch (20) was used to subdivide the grouping, with a cutoff of 80% nucleotide identity, to account for randomly overlapping UMIs. Each of the resulting groupings is assumed to represent reads arising from a single RNA. Reads within each grouping were then aligned, and a consensus sequence determined. Finally, duplicate reads were collapsed into a single processed sequence for analysis. Collapsing duplicate reads ensures that each processed sequence represents a sequence from a single B cell and our analysis is not confounded by expression level.


Input files:

1. The read 
2. The primers sequences available online at the table below.


Output files:

1. sample_collapse-unique.fastq
2. log and log tab file for each step.

Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

1. Assemble pairs
2. FilterSeq quality
3. MaskPrime - Extracting primers and UMIs
4. ClusterSets - Cluster sequences based on sample barcode, UMI, and constant region primer
5. BuildingConsensus
6. CollapseSeq - Collapsing duplicates reads


Primers used:

* [BIOMED2_SPRIMERS](https://github.com/PipeAIRR/pipeline4/blob/master/primers/S_primer.fasta)

* [BIOMED2_CPRIMERS2](https://github.com/PipeAIRR/pipeline4/blob/master/primers/C_primers_file.fasta)

* [BIOMED2_VHFR1](https://github.com/PipeAIRR/pipeline4/blob/master/primers/V_primers_file.fasta)
