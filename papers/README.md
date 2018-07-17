# Hi-C papers

The latest on top

## Methodological Reviews

- Nicoletti, Chiara, Mattia Forcato, and Silvio Bicciato. “Computational Methods for Analyzing Genome-Wide Chromosome Conformation Capture Data.” Current Opinion in Biotechnology 54 (December 2018): 98–105. https://doi.org/10.1016/j.copbio.2018.01.023. - 3C-Hi-C tools review, Table 1 lists categorizes main tools, Figure 1 displays all steps in technology and analysis (alignment, resolution, normalization, including accounting for CNVs, A/B compartments, TAD detection, visualization). Concise description of all tools.

- Waldispühl, Jérôme, Eric Zhang, Alexander Butyaev, Elena Nazarova, and Yan Cyr. “Storage, Visualization, and Navigation of 3D Genomics Data.” Methods, May 2018. https://doi.org/10.1016/j.ymeth.2018.05.008. - Review of tools for visualization of 3C-Hi-C data, challenges, analysis (Table 1). Data formats (hic, cool, BUTLR, ccmap). Database to quickly access 3D data. Details of each visualization tool in Section 4

- Forcato, Mattia, Chiara Nicoletti, Koustav Pal, Carmen Maria Livi, Francesco Ferrari, and Silvio Bicciato. “Comparison of Computational Methods for Hi-C Data Analysis.” Nature Methods, June 12, 2017. https://doi.org/10.1038/nmeth.4325. - Hi-C processing and TAD calling tools benchmarking, [Table 1](https://www.nature.com/articles/nmeth.4325/tables/1), simulated (Lun and Smyth method) and real data. Notes about pluses and minuses of each tool. TAD reproducibility is higher than chromatin interactions, increases with larger number of reads. Consistent enrichment of TAD boundaries in CTCF, irrespectively of TAD caller. Hi-C replication is poor, just a bit more than random. Supplementary table 2 - technical details about each program, Supplementary Note 1 - Hi-C preprocessing tools, Supplementary Note 2 - TAD callers. Supplementary note 3 - how to simulate Hi-C data. Supplementary note 6 - how to install tools. https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf

- Yardımcı, Galip Gürkan, and William Stafford Noble. “Software Tools for Visualizing Hi-C Data.” Genome Biology 18, no. 1 (December 2017). https://doi.org/10.1186/s13059-017-1161-y. - Hi-C technology, data, and visualization review. Suggestion about graph representation.

- Ay, Ferhat, and William S. Noble. “Analysis Methods for Studying the 3D Architecture of the Genome.” Genome Biology 16 (September 2, 2015): 183. https://doi.org/10.1186/s13059-015-0745-7. - Hi-C technology and methods review. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7#Tab1) - list of tools. Biases, normalization, matrix balancing. Extracting significant contacts, obs/exp ratio, parametric (powerlaw, neg binomial, double exponential), non-parametric (splines). 3D enrichment. References. TAD identification, directionality index. Outlook, importance of comparative analysis

## General Reviews

- Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation.” Molecular Systems Biology 11, no. 12 (December 23, 2015): 852–852. doi:10.15252/msb.20156492. http://msb.embopress.org/content/msb/11/12/852.full.pdf - 3D genome organization parts. Well-written and detailed. References. Technologies: FISH, 3C. 4C, 5C, Hi-C, GCC, TCC, ChIA-PET. Typical resolution - 40bp to 1Mb. LADs - conserved, but some are cell type-specific. Chromosome territories. Cell type-specific. inter-chromosomal interactions may be important to define cell-specific interactions. A/B compartments identified by PCA. Chromatin loops, marked by CTCF and Cohesin binding, sometimes, with Mediator. Transcription factories

- Dekker, Job, Marc A. Marti-Renom, and Leonid A. Mirny. “Exploring the Three-Dimensional Organization of Genomes: Interpreting Chromatin Interaction Data.” Nature Reviews. Genetics 14, no. 6 (June 2013): 390–403. https://doi.org/10.1038/nrg3454. https://www.nature.com/articles/nrg3454 - 3D genome review. Chromosomal territories, transcription factories. Details of each 3C technology. Exponential decay of interaction frequencies. Box 2: A/B compartments (several Mb), TAD definition, size (hundreds of kb). TADs are largely stable, A/B compartments are tissue-specific. Adjacent TADs are not necessarily of opposing signs, may jointly form A/B compartments. Genes co-expression, enhancer-promoters interactions are confined to TADs. 3D modeling.

## Comparative analysis of Hi-C data

- Bonev, Boyan, Netta Mendelson Cohen, Quentin Szabo, Lauriane Fritsch, Giorgio L. Papadopoulos, Yaniv Lubling, Xiaole Xu, et al. “Multiscale 3D Genome Rewiring during Mouse Neural Development.” Cell 171, no. 3 (October 2017): 557–572.e24. https://doi.org/10.1016/j.cell.2017.09.043. - Mouse development 3D genomics. Ultra-deep Hi-C during mouse neural differentiation. Transcription correlated with 3D structure, but insufficient to change it. During differentiation, active TADs increase interaction strength within them, inactive - decrease. Polycomb binding is disrupted with differentiation. Number of TAD boundaries decreases with differentiation, correspondingly, the size of TADs increases. More boundaries without CTCF binding. Intro into A/B compartments, TADs. Biological replicates data, highly correlated at all resolutions. Used Shaman for randomization, reference to (Olivares-Chauvet et al., 2016) Nature paper for methods https://www.nature.com/nature/journal/v540/n7632/full/nature20158.html.

- Dixon, Jesse R., Inkyung Jung, Siddarth Selvaraj, Yin Shen, Jessica E. Antosiewicz-Bourget, Ah Young Lee, Zhen Ye, et al. “Chromatin Architecture Reorganization during Stem Cell Differentiation.” Nature 518, no. 7539 (February 19, 2015): 331–36. https://doi.org/10.1038/nature14222. - Reorganization of the 3D genome during human embryonic stem cell differentiation. ~36% of active/inactive A/B compartments change. Fig. 1 - how to compare PC1 vectors. Gene expression changes in compartments that switch states (up/down in B to A vs. A to B swiching), small but significant. Changes in _interaction frequency_ (increase/decrease - switch from B to A vs. A to B), correlated positively with activating epigenomic marks H3K27ac, DNAse hypersensitive sites, CTCF, negatively with repressive marks H3K27me3 and H3K9me3. Random Forest to select epigenomic features classifying changes in interaction frequency - H3K4me1 enhancer mark. Haplotype-resolved Hi-C maps are similar. Methods: A/B compartment visualization, assignment to active/repressed states by gene density. Overlap-based comparison, but p-value for "concerted" changes in IFs within one domain.




## Thesis

- "4D Nucleome of Cancer" by Laura Seaman, 2017, https://deepblue.lib.umich.edu/handle/2027.42/140814

- "GENOME ANALYSIS IN THREE DIMENSIONS: FUNCTIONAL ANALYSIS OF HI-C DERIVED DATASETS" by Robert Sugar, 2014, https://www.ebi.ac.uk/sites/ebi.ac.uk/files/shared/documents/phdtheses/RobertThesis_2014-12-12_v03_CORRECTED.pdf. 3D overview, TADs, CTCF, Cohesin, Mediator. 3D, epigenomics and transcription regulation. Technical aspects of 3C technologies, alignment. Existing pipelines (Hicpipe, hiclib, Hi-Five recimmended). Methods for simulating Hi-C data. Some methods are implemented in GOTHiC R package. Capture Hi-C technology for targeted 3C capturing. Real-life experiments, analysis, detection of significant interactions, enrichment in epigenetic marks, association of chromatin interactions with gene expression, GWAS SNP enrichment in the promoters of genes. Future of 3C technologies.

## Misc

- A list of tools available for data analysis and/or visualization of 4DN-related datasets. https://www.4dnucleome.org/software.html

- `manual_180319.docx` - "Genome-Assembly-Cookbook", https://github.com/theaidenlab/Genome-Assembly-Cookbook