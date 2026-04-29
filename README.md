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
