# 142 Y chromosomes
Code repo pertaining to a subset of analyses of the Pan-Y (joint HGSVC + HPRC) project (ignore file numbering as things change during analyses).

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

## Constraint:
* Yq12_Constraint_Testing_Submission.ipynb - Tests the evolutionary dynamics of the two major Yq12 repeat families, DYZ1 and DYZ2, using phylogenetically informed models. The notebook evaluates whether changes in DYZ1 and DYZ2 are evolutionarily coupled, whether their evolutionary rates differ, and whether these patterns are supported by the Y-chromosome phylogeny.
* AZFc_Constraint_Testing_Submission.ipynb - Tests whether the observed diversity of human AZFc architectures is more structurally constrained than expected under a permissive model of NAHR-mediated rearrangement. The notebook reconstructs architecture changes on the Y phylogeny and uses continuous-time simulations of inversions, deletions, and duplications to compare the observed range of AZFc structures with the null expectation.
* Centromere_Constraint_Testing_Submission.ipynb - Tests whether the functional Y centromere domain is constrained relative to variation in the surrounding DYZ3 α-satellite array. The notebook evaluates phylogenetic signal, scaling between CDR and total DYZ3 array size, and whether CDR size is less variable than expected under evolutionary models fitted to the surrounding satellite array.

## DAZ_Submission:
* DAZ_PSV_IDWork-part2.ipynb — Parses DAZ paralogous sequence variants, maps PSV alleles to DAZ1–DAZ4 using BLAT hits, and summarizes/visualizes PSV-based DAZ gene identity patterns.
* DAZ_PSVs_patterns_part3.ipynb — Uses PSV placements and RepeatMasker context to identify recurring DAZ exon/RRM motif patterns and compare pre- versus post-LINE sequence structure.
* DAZGenes_KMERs_part4.ipynb — Mines DAZ gene sequences for unique and multi-hit k-mers, then visualizes shared k-mer relationships among DAZ genes across samples.
* DAZGenes_KMERs-wNTVisualization_part5.ipynb — Extends the DAZ k-mer analysis by mapping shared k-mers back onto nucleotide positions and visualizing their locations across paired DAZ gene sequences.
* Part1_Publication_DAZ_Exons_AllSamples-HPRC-HGSVC3-CEPH_DAZNAMES-wRepeatMaskerInformation-additionalExonInformatoin-Submission.ipynb — Builds publication-ready DAZ exon annotations across HPRC/HGSVC/CEPH samples, integrating DAZ names, exon structure, RepeatMasker information, haplogroups, and phylogenetic ordering for visualization.

## DNMs_Yq12:
* Visualize_DNM_Yq12_SNVs-Submission.ipynb — Integrates validated Yq12 de novo SNVs with repeat architecture, donor metrics, and father–son pair metadata to generate publication-ready visualizations of SNV positions, clustering, and candidate gene conversion events across Yq12.
* DNM_GeneConversion_DonorCheck_Submission_04232026.ipynb — Implements the full seed-and-extend donor-discovery and statistical classification pipeline for Yq12 de novo SNVs, identifying homologous donor loci, testing gene conversion support across homology thresholds, collapsing concordant SNV clusters, and classifying events as likely de novo mutation or gene conversion.

## DYZ17_DYZ19_Yq12:
* Repeat_Array_Lengths.ipynb — Quantifies and compares lengths of major Y-chromosome repeat arrays (including Yq12, centromeric, and ampliconic structures) across assemblies, summarizes haplogroup-level size variation, and generates comparative publication-quality plots of repeat architecture diversity.

## GeneAnnotation:
* PullAmpliconicGenesSequences.ipynb — Extracts and organizes full ampliconic gene family sequences (e.g., DAZ, RBMY, TSPY, BPY2) from Y-chromosome assemblies for downstream comparative annotation and copy-number analyses.
* Part2_VisualizeAmpliconicGeneCopyNumbers.ipynb — Aggregates ampliconic gene annotations across samples, quantifies gene-family copy numbers, and visualizes structural variation in major multicopy Y-linked genes across haplogroups.
* Part3_Combine_All_Gene_Annotations-AllContigs.ipynb — Merges gene annotations from all Y contigs and samples into a unified dataset for cross-sample gene presence, structure, and positional analyses.
* Part3_Combine_All_Gene_Annotations-BEDFILE_allContigs(2).ipynb — Converts merged multi-contig gene annotations into standardized BED-style interval datasets for downstream genomic analyses and visualization.
* Part4_non_Y_Genes.ipynb — Identifies and catalogs non-Y gene annotations or off-target gene mappings present within assemblies to separate true Y-linked genes from potential annotation artifacts.
* part5_ampliconic_pseudogenes.ipynb — Detects, classifies, and visualizes pseudogenized ampliconic gene copies, helping distinguish intact versus disrupted multicopy gene family members across Y assemblies.

## MEI_Submission:
* Part5_Phylogeny.ipynb — Constructs and visualizes Y-chromosome phylogenetic relationships across assemblies using structural and sequence-based variation, generating ordered haplogroup frameworks for downstream comparative analyses.
* FindAllTEs_Elements_Part1.ipynb — Identifies transposable element insertions across Y assemblies by parsing repeat annotations, extracting candidate TE loci, and organizing element-specific insertion calls.
* FindAllTEs_part2.ipynb — Refines TE insertion discovery by consolidating, comparing, and classifying transposable element calls across samples to identify shared and lineage-specific insertions.
* part4_FilterCalls_WholeY.ipynb — Applies quality-control, positional, and repeat-context filters to whole-Y TE or structural variant callsets to generate a high-confidence final dataset for evolutionary analysis.

## RBMY:
* Part1_PullRBMYSequences-SELECTED_GOOD-PUBLICATIONWORTHYexactMatches.ipynb — Extracts high-confidence RBMY gene sequences from assemblies, identifies exact-match RBMY copies, and prepares curated sequence sets for structural and evolutionary analyses.
* Part2_ReadInAlignments-NetworkAnalaysis-visualizeRBMY-WorkingPlainOriginal.ipynb — Reads RBMY multiple-sequence alignments, quantifies sequence similarity, builds network-based RBMY relationship models, and visualizes structural organization across RBMY copies.
* Part2_ReadInAlignments-visualizeRBMY-withUpstreamClusters.ipynb — Extends RBMY alignment visualization by incorporating upstream regulatory or flanking sequence clusters to better define proximal/distal RBMY structural groupings.
* Part3_TFMotif_RBMY_Locations.ipynb — Identifies transcription factor motif locations near RBMY genes, integrates promoter/regulatory annotations, and visualizes motif distributions relative to RBMY structural organization.
* RBMY_Groups_Numbers.ipynb — Quantifies RBMY copy-number variation and structural grouping patterns across samples, summarizing proximal/distal array organization and haplogroup-level diversity.

## TSPY:
* Part1_PullTSPYSequences-TF_Submission.ipynb — Extracts curated TSPY gene and upstream regulatory sequences from Y assemblies for publication-quality comparative analyses of TSPY structural diversity.
* Part2_ReadInAlignments-visualize-Networks-UpstreamMotifs_Submission.ipynb — Reads TSPY sequence alignments, constructs network relationships among TSPY copies, integrates upstream motif architecture, and visualizes array organization across samples.
* Part2_ReadInAlignments-visualize-pseudogenes_OnlyForVisualization.ipynb — Visualizes pseudogenized TSPY copies and their sequence divergence patterns to distinguish intact versus disrupted TSPY family members.
* Part2_ReadInAlignments-visualize-Submission_OnlyForVisualization.ipynb — Produces publication-ready visualizations of TSPY sequence alignments, structural clusters, and tandem array organization across diverse Y assemblies.
* Part7_TFMotif_LOcations.ipynb — Maps transcription factor binding motifs across TSPY loci and surrounding regulatory regions, integrating motif architecture with TSPY structural variation.

## Y_Structure:
* Ychromosome_Comparison.ipynb — Performs large-scale comparative analyses of complete Y-chromosome assemblies by integrating structural annotations, repeat architecture, gene content, and haplogroup diversity to identify conserved versus variable genomic features across human populations.

## Yq12_Submission:
* Yq12_Structure_and_Inversions.ipynb — Characterizes large-scale Yq12 repeat architecture across assemblies, identifies structural configurations, and analyzes inversion patterns within DYZ1/DYZ2 heterochromatic arrays.
* Yq12_InversionBreaks-AllSamples-ProximalInversion.ipynb — Detects and compares proximal Yq12 inversion breakpoints across all samples, focusing on breakpoint consistency, orientation, and repeat-context structure.
* Yq12_InversionBreaks-AllSamples.ipynb — Systematically identifies Yq12 inversion breakpoints across the full cohort, compares breakpoint architectures, and classifies recurrent versus lineage-specific inversion events.
* Yq12_InversionBreaks-BreakPoints_Submission.ipynb — Produces publication-ready analyses of Yq12 inversion breakpoint intervals, integrating sequence alignments, repeat annotations, and structural classifications.
* Yq12_InversionBreaks-Copy1.ipynb — Earlier developmental version of the Yq12 inversion breakpoint pipeline used to identify candidate breakpoint regions and structural rearrangement signatures.
* Yq12_InversionBreaks.ipynb — Core workflow for detecting Yq12 inversion breakpoints through sequence comparison, breakpoint-window discovery, and repeat-structure analysis.

## hmmerFilterCode:
* filterYq12HMMER.py — Filters and consolidates nhmmer hits for major Yq12 repeat families (DYZ1, DYZ2, DYZ18, and associated heterochromatic elements), merges adjacent repeat blocks, resolves overlapping annotations, trims conflicting regions, and produces refined Yq12 structural architecture coordinates for downstream comparative analyses.
* filterDYZ19HMMER.py — Processes DYZ19/centromeric nhmmer outputs to identify, merge, and refine high-confidence DYZ19 array blocks, estimate total repeat copy numbers, detect structural disruptions such as possible inversions, and summarize centromeric array organization per assembly.
* filterTSPYHMMER.py — Parses TSPY nhmmer outputs to reconstruct full TSPY gene arrays, merge fragmented or overlapping TSPY hits, distinguish major tandem arrays from isolated copies such as TSPY2, identify pseudogene or fragmented structures, and generate refined TSPY architectural annotations across assemblies.
