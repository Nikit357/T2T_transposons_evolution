
Decision Letter (MBE-25-1299)

From:

EiC.MBE@gmail.com, eassist.mbe@publishingsolutions.net

To:

danya.nikitin.orel@gmail.com

CC:

EiC.MBE@gmail.com, eassist.mbe@publishingsolutions.net

Subject:

Editorial Decision to Reject MBE-25-1299

Body:

02-Feb-2026

MS: MBE-25-1299
Title: Transposable element-host genome evolutionary arms race revealed by multi-modal epigenomic profiling in a telomere-to-telomere human genome reference

Dear Dr. Nikitin:

Thank you for submitting your manuscript to Molecular Biology and Evolution (MBE). We regret to inform you that it did not receive high enough priority for publication after an in-depth review by the editors and the peer reviewers. Specific comments from the editors and/or external reviewers are included below.

In general, MBE seeks to publish research, methods, and resources of broad significance in molecular evolutionary biology. Even when the external reviewers find a manuscript to be scientifically and technically sound, the ultimate priority for publication is determined based on the novelty and impact of the work presented. MBE does not publish manuscripts judged by the reviewers to contain mostly descriptive work, confirmatory results, and discoveries with a limited gene and taxonomic scope.  All of these factors were considered in deciding the publication priority for your manuscript.

Thank you for considering Molecular Biology and Evolution, and please continue to consider MBE as a venue for the publication of your best work.

Sincerely,

Board of Editors
Molecular Biology and Evolution

Associate Editor
Editors’ comments to the author:
The comments and recommendations from two expert reviewers are now available for your manuscript. The reviewers assessed the reported discoveries as having medium significance and medium potential scientific impact. Both reviewers provided thoughtful suggestions and identified critical technical issues that need to be thoroughly addressed. Specifically, the estimation of enrichment of repressive marks at repetitive TEs and the estimation of TE age based on RepeatMasker tracks both require complete revision to ensure the robustness of the results. Given the substantial analytical work required to address these concerns, significant new analyses will be necessary to bring the manuscript to a standard acceptable for publication in MBE.

Reviewer(s)' comments to author:
Reviewer: 1

Comments to the Author
In this manuscript, the author investigates the association between the age of transposable elements and their epigenetic marks, using the telomere-to-telomere (T2T) human genome assembly and ENCODE dataset. The goals are to study the evolutionary and epigenetic dynamics of different TEs and evaluate the relative contributions of different epigenomic marks on TE dynamics. The question is interesting and significant. But I am uncertain about the robustness of the results and how this work is distinguished from previous works cited.

Major comments
1. To estimate epigenetic enrichment around TEs, the author used the enrichment tracks computed from the T2T ENCODE ChIP-seq pipeline. The average enrichment value of the genomic interval for each TE was calculated. Because of the repetitiveness of TE sequences, the mappability of the ChIP-seq reads in the TE regions could affect the estimation. This is more challenging for recently active TE families with many similar copies in the genome, increasing the measurement errors. I feel it would be important to check the mappability for different TE families, especially for those showing age-dependent epigenetic marks changes, like SVA and Alu.
2. I am not sure that the statistical methods are well-justified. First, while Pearson correlation coefficient can be informative for the magnitude of dependencies between age and epigenetic activity, the CIs and p-values (Table 1 & 2) are not reliable if the normality and homoscedasticity assumptions are violated. SVA and SINE seem to have violated the assumptions. I would also suggest trying Spearman’s correlation test and using the coefficient rho for subsequent analysis. Second, the author found a negative relation between the Pearson correlation coefficient and copy number among families and fit a sigmoid model (Fig 3A). To account for low-copy families inflated by random error, a one-time random permutation between divergence and enrichment was performed (line 272-274). Why only use one-time? Moreover, F-test was used to compare the observed sigmoidal fit and the fit on permutated data. I find it hard to interpret because it suggests a sigmoid model fits the observed data better than random; it doesn’t tell whether the relation is more negative than random. Perhaps a better way to test is to calculate Spearman correlation coefficient instead of a sigmoid fit and test whether the coefficient is more negative than a random expectation based on many permutations.

Other comments
1. In the method section, line 633-635, the author mentioned briefly that the TE annotations were obtained from the RepeatMasker track. Perhaps more introduction for the source of annotation is helpful? Is it from the “Hoyt et al. 2022, From telomere to telomere: The transcriptional and epigenetic state of human repeat elements” study?
2. One of the major advantages of using the T2T assembly for TE analysis is that it improves the annotations of active TEs in heterochromatin regions. Maybe it’s worth exploring whether the patterns differ for TEs located in pericentromeric regions? How about TEs located in genes?
3. Fig 1D indicates that for divergence scores higher than certain values (~220), it is not a good indicator for age inferred with a phylogenetic approach. Perhaps filter out TEs with divergence score above certain values (e.g., 200, 220, and 250) and see whether the results are consistent?
4. Which specific results were generated using RetroSpect? Could you briefly explain how RetroSpect works?
5. Is SVA a TE family belonging to the SINE class? Fig 2 seems to imply that SVA is a class different from SINE. Also, does the DNA class include Helitrons? If not, maybe use TIR instead of DNA?
6. Line 314, “All TEs classes exhibited statistically significant differences” means different from what? Did you perform pairwise comparisons?

Reviewer: 2

Comments to the Author
I have now reviewed the manuscript of Daniil Nikitin entitled "Transposable element–host genome evolutionary arms race revealed by multi-modal epigenomic profiling in a telomere-to-telomere human genome reference". In this work, the author present an analysis of the correlation between ChIP-seq enrichment signal for 7 epigenetic marks (CTCF, H3K4me1, H3K4me3, H3K9ac, H3K27ac, H3K27me3, H3K9me3) and transposable element (TE) relative evolutionary age, approximated by their divergence from their consensus sequence, as reported by RepeatMasker on the T2T reference genome CHM13.

As presented by the author in introduction, this work provide an opportunity to revisit the hypothesis that the association between TE and epigenetic regulation is dynamic through time and can evolve in the context of a continuous arms race between TE (which success can be measured in the ability to produce new copies) and host genome genome defense mechanisms which evolve to limit repress potentially damaging TE activity.

I generally concur with the author regarding the cautious interpretation of the results. My principal reservation is however related to the direct use of the divergence value provided by the Repeatmasker track (as detailed below), but I provide a practical solution to it. Other comments and suggestions that I believe can strengthen the manuscript are also provided. Perhaps, overall I would refrain to extrapolate the results as evidence of arms race but rather as supporting the hypothesis of ongoing arms race.

1. Divergence

The author cite Giordano et al. (2007) to indicate that "The average divergence from the consensus sequence of each TE family was used as a proxy for evolutionary age", however it seem to me that the point of the cited publication is to indicate that divergence, as reported in the Repeatmasker table is an imperfect metric. In their discussion, the author states: "Because the rate of sequence divergence (the molecular clock) may not be constant over time or between lineages, the age estimates of TEs based on percent divergence may not be entirely reliable, especially for the older, more diverged elements. Our method to determine relative ages of TEs is not dependent on the percent divergence from derived consensus sequences or on an assumption of a constant molecular clock, and hence can be applied to all TEs in a given genome that have interacted with (inserted into or been interrupted by) enough TEs."

In order to validate the divergence as good proxy for the analysis, the author use data from Kosuge, Ito and Hamada (2024) to test the correlation between divergence and TE age estimation. Accordingly, I think the statement "Collectively, these results support the use of average divergence as a reliable proxy for the evolutionary age of TE families and subfamilies" might be revised to indicated that this appear true up to 110 Myo with the current data.

RepeatMasker however is able to provide a somewhat corrected estimate of the molecular divergence through the Kimura 2-parameters (K2P) value, which is available in the .align outputs of standard runs. I am aware though that these data were not directly available through the UCSC genome browser, and I have thus arranged to make them available by the RepeatMasker team who worked on the T2T paper (Nurk et al., 2022). The author can now find it at: https://www.repeatmasker.org/genomes/hs1/rmsk4.2.2_dfam3.9_rb20181026_rmb/hs1.fa.align.gz -- this .align file correspond to the .out file similar to the table used in this manuscript (also available here: https://www.repeatmasker.org/genomes/hs1/rmsk4.2.2_dfam3.9_rb20181026_rmb/hs1.fa.out.gz)

For each hit, the .align file shows a line similar to the .out, followed by the alignment followed by a few statistics including the line "Kimura (with divCpGMod) = XXX". To obtain the K2P divergence value in the .align file corresponding to the copy in the .out file, one can use the hit ID value, which is given in the last column of the .out file (and reproduced at the end of the hit line in the .align file -- see also: https://www.repeatmasker.org/webrepeatmaskerhelp.html). Note that several lines in the .out can share the same hit ID if a copy is interrupted by a more recent element. Multiple alignments segments will also exist in the .align. In this case, the divergence/K2P value can be calculated by computing the weighted average, using the length of each sub hit (each segment). Finally, be aware of a special case regarding LINE element. Since RepeatMasker switched to DFAM, LINE elements in the human lineages have been split in 2 to 3 segments (5' end, 3' end and orf2) in order to improve the precision of RepeatMasker however, this is not shown in the .out table. So, for a given hit ID corresponding to a LINE element, you will likely find 2 or more alignments matching the 3'end and orf2 (perhaps 5' end if full-length hit) in the .align file. It is also possible that the LINE segments in the .align do not match the final family or subfamily reported in the .out: usually, it is the name of the 3'end fragment that will be reported, as it carries diagnostic SNP for the subfamilies. If you have questions about how to interpret these files, please reach out to help@dfam.org


2. I suggest to add some details about how the ENCODE enrichment values have been calculated for context, as they are underlying the whole analysis.

3. Please also add definitions of Class, family and subfamily as used in the context of the paper to help non-specialist of human TE.

4. Overall, I concur with the statement of the author that the data "suggest[ing] that individual families within these numerically dominant classes may nevertheless harbor similarly pronounced regulatory dynamics that are masked at the aggregate class level." It seem to me that class-level correlation are observed in classes that have bi-modal distributions. I hope that using K2P can improve these correlations.

5. SINE and SVA clustering in Figure 2: could that be explained by homology of the SINE domain in the SVA with Alu? Perhaps they are too different.

6. Family-level correlations. I have a few reservation with the sigmoid fit with absolute Pearson's correlation coefficients and copy numbers:
- Can you detail how the TE copy number were counted: were they tabulated by counting the lines or hit IDs?
- It seems that the fit is mainly driven by DNA families, that are old and may be in fact degenerate. Perhaps confusing low copy number with extreme age.

7. In general Class, families of subfamilies showing best correlations with a mark for a given cell type have a bi-modal divergence. I would like to see, for at least some of the showcased examples in figure 2, 3 and 4, if there are significant difference in copy length between the youngest (leftmost) and oldest (rightmost) modes as well as pileups of the genomic hits over the family/subfamily consensus (figures 3 and 4) for the younger and older copies. The reasoning is that the epigenetic marks are not necessarily deposited homogeneously along the TE sequence, and older TE instances (with high divergence value) may reflect shorter fragments that differentially include or not binding sites rather than having lost or gained binding side through evolution. Say an older SVA_B fragment is seen with less H3K9me3 levels than a young SVA_B instance, but perhaps the old SVA_B instance represent a fragment of the TE that even young instance do not get typically deposited with H3K9me3.

Date Sent:

02-Feb-2026
