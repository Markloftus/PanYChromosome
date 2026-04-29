# 142 Y chromosomes
Code repo pertaining to a subset of analyses of the Pan-Y (joint HGSVC + HPRC) project.

## AZFc_Submission:
* ColorBlock_Part2.5_06112025_compareArang.ipynb — Compares refined AZFc color-block annotations against BED annotations (produced by Arang) and appends matching annotation information to each sample’s color-block dataframe.
* FindBreakpoints-Part1_Minimap2_ColorBlocks_pullSequences_BREAKPOINTS.ipynb — Builds color-block groupings and pulls FASTA sequences for paired samples so corresponding AZFc blocks can be aligned for breakpoint discovery.
* FindBreakpoints-PartA_b2b3_Inversion.ipynb — Extracts front/back sequences around candidate b2/b3 inversion breakpoints for selected sample pairs.
* FindBreakpoints-PartA_DELETIONS.ipynb — Uses shared unique k-mers between deleted and full-length AZFc haplotypes to visually compare candidate deletion breakpoint regions.
* FindBreakpoints-PartA_T2T-vs-hg38_Inversion.ipynb — Extracts candidate breakpoint-flanking sequences for T2T-like versus GRCh38-like AZFc inversion comparisons.
* FindBreakpoints-PartB_DELETIONS.ipynb — Pulls deleted-sample and full-haplotype flanking sequences into FASTA files for deletion breakpoint alignment.
* FindBreakpoints-PartB_Read_in_T2TvsGRCh38-Alignments.ipynb — Reads multiple-sequence alignments from T2T-vs-GRCh38 inversion breakpoint FASTAs and plots front/back breakpoint match patterns.
* FindBreakpoints-PartC_DELETIONS_Read_Alignments.ipynb — Reads deletion breakpoint alignments, maps alignment columns back to genomic coordinates, and visualizes deleted versus full haplotype breakpoint-supporting matches.
* FindBreakpoints-PartC_Read_in_b2b3Inversion.ipynb — Reads b2/b3 inversion breakpoint alignments, identifies breakpoint windows, and annotates them with repeat information.
* FindBreakpoints-PartC_Read_in_grrgInversion.ipynb — Reads gr/gr inversion breakpoint alignments, identifies candidate breakpoint windows, and summarizes repeat/block context for those windows.
* FindBreakpoints-PartD_DELETIONS_FindWindows.ipynb — Processes deletion alignment details to define candidate breakpoint windows between conserved alignment blocks.
* FindBreakpoints-PartD_INVERSIONS.ipynb — Groups inversion breakpoint windows by color block, orientation, and relative overlap, then tests repeat enrichment around inversion breakpoints.
* FindBreakpoints-PartE_DELETIONS.ipynb — Groups deletion breakpoint windows by relative position/color block and tests repeat-class enrichment within deletion breakpoint intervals.
* FindBreakpoints-PartF_DELETIONS.ipynb — Further collapses and compares deletion breakpoint windows by color block and orientation to identify shared breakpoint-region groups.
* Pullhg38_ColorRegion_Part2_09102025.ipynb — Uses GRCh38 color-region k-mers to refine color-block coordinates across sample assemblies.
* Pullhg38_ColorRegion_windows-Part1_AllSamplesRun-09102925-original.ipynb — Generates initial GRCh38-derived color-region windows across all samples using k-mer density and interval merging.
* VisualizeColorBlock_overview_Part4.ipynb — Creates overview visualizations of AZFc color-block architecture across samples, integrating QC regions, gene annotations, copy-number data, and phylogenetic ordering.

## DAZ_Submission:
* DAZ_PSV_IDWork-part2.ipynb — Parses DAZ paralogous sequence variants, maps PSV alleles to DAZ1–DAZ4 using BLAT hits, and summarizes/visualizes PSV-based DAZ gene identity patterns.
* DAZ_PSVs_patterns_part3.ipynb — Uses PSV placements and RepeatMasker context to identify recurring DAZ exon/RRM motif patterns and compare pre- versus post-LINE sequence structure.
* DAZGenes_KMERs_part4.ipynb — Mines DAZ gene sequences for unique and multi-hit k-mers, then visualizes shared k-mer relationships among DAZ genes across samples.
* DAZGenes_KMERs-wNTVisualization_part5.ipynb — Extends the DAZ k-mer analysis by mapping shared k-mers back onto nucleotide positions and visualizing their locations across paired DAZ gene sequences.
* Part1_Publication_DAZ_Exons_AllSamples-HPRC-HGSVC3-CEPH_DAZNAMES-wRepeatMaskerInformation-additionalExonInformatoin-Submission.ipynb — Builds publication-ready DAZ exon annotations across HPRC/HGSVC/CEPH samples, integrating DAZ names, exon structure, RepeatMasker information, haplogroups, and phylogenetic ordering for visualization.

## DNMs_Yq12
* Visualize_DNM_Yq12_SNVs-Submission.ipynb — Integrates validated Yq12 de novo SNVs with repeat architecture, donor metrics, and father–son pair metadata to generate publication-ready visualizations of SNV positions, clustering, and candidate gene conversion events across Yq12.
* DNM_GeneConversion_DonorCheck_Submission_04232026.ipynb — Implements the full seed-and-extend donor-discovery and statistical classification pipeline for Yq12 de novo SNVs, identifying homologous donor loci, testing gene conversion support across homology thresholds, collapsing concordant SNV clusters, and classifying events as likely de novo mutation or gene conversion.
