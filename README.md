# Vitamaize_NAM_GWAS_LabVersion
Contains all scripts for tocochromanol project, proofed August 2016 and October 2017.
All scripts written by C.H. Diepenbrock and C.B. Kandianis, with contributions from A.E. Lipka, unless otherwise specified. Please reference Materials and Methods for more details. Scripts are sectioned as follows:

(1) Joint Linkage (JL) Analysis
Core Analysis: determine trait-specific significance thresholds through permutations and conduct JL. These analyses used TASSEL (Bradbury et al. 2007) with minor adaptation.
1a - check for multicollinearity. Modified TASSEL 3 source build was then used within NetBeans to fit a JL model using user-specified markers (those remaining after multicollinearity correction) for generation of final model and support intervals.
1b - generate JL model effect estimates (from transformed BLUEs) and residuals for GWAS.
1c - generate JL model effect estimates (from untransformed BLUEs) to discern magnitude and direction of QTL effects.

(2) Genome-Wide Association Study (GWAS)
Core Analysis: determine trait-specific significance thresholds through permutations and conduct GWAS. These analyses were scripted with minor adaptation from Dataset S3 in Wallace et al. (2013), and use TASSEL 4.
2a - merge GWAS results from each chromosome into a single file per trait, and (optionally) visualize by cM window.

(3) JL-GWAS Summaries
3a - calculate PVE, compile tabular summary.
3b - merge physically overlapping individual-trait support intervals into common support intervals.
3c,d,e - intersect JL and GWAS results, for tabular summary and in preparation for triangulation.
3f - calculate LOD scores.

(4) Triangulation: calculating correlations 
4a,b,c - extract GWAS variant genotype state scores. Script 4b labeled 'INDIV_QTL' is to be used when genotype input files exceed 1 Mb.
4d and e - calculate pairwise correlations and P-values between RNA-seq expression levels, JL effect estimates (from transformed BLUEs), and GWAS variant genotype state scores within each individual-trait support interval.
4f and g - compile triangulation correlation output at the level of each common support interval.

(5) Pleiotropy
5a - calculate correlations of effect estimates when applying the final JL model for each trait to every other trait, with significance testing. This analysis was as conducted in Brown et al. (2011) with minor adaptation.
5b - re-format 5a results for visualization and in preparation for plotting.
5c - plot trait network depicting magnitude and direction of shared genetic basis across traits. An additional minimum level of proportion of shared QTL required for plotting can be imposed at this stage.

(6) Epistasis
6a - Permutate phenotypes by family to establish significance thresholds; test additive x additive terms; calculate PVE for terms showing significance.
6b - generate JL model effect estimates (from untransformed BLUEs).

(7) Master gene summaries
Plotting scripts to generate Figures 2 and S5, which show various aspects of the triangulation results for the 14 NAM JL-QTL in which genes were identified. These scripts were used as written by Dan Ilut with minor adaptation.

References:
Bradbury et al. 2007. "TASSEL: software for association mapping of complex traits in diverse samples." Bioinformatics 23(19):2633-2635. doi: 10.1093/bioinformatics/btm308.
Brown et al. 2011. "Distinct Genetic Architectures for Male and Female Inflorescence Traits of Maize." PLOS Genetics 7(11): e1002383. doi:10.1371/journal.pgen.1002383.
Wallace et al. 2013. "Association Mapping across Numerous Traits Reveals Patterns of Functional Variation in Maize." PLOS Genetics 10(12): e1004845. doi:10.1371/journal.pgen.1004845.
