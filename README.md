# A collection of tools and papers related to Hi-C data analysis

Slowly growing as notes from my Zotero collection are getting organized. A related repository holds references to Hi-C data, https://github.com/mdozmorov/HiC_data. Issues and/or Pull requests to add other data are welcome!

# Table of content

* [Pipelines for Hi-C data processing](#pipelines)
* [Normalization of Hi-C data](#normalization)
* [Reproducibility and QC of Hi-C data](#reproducibility)
* [Significant interaction (peak) callers](#interaction-callers)
* [TAD callers](#tad-callers)
* [Prediction of 3D features](#prediction-of-3d-features)
* [SNP-oriented Hi-C analysis](#snp-oriented)
* [Structural variant detection](#structural-variant-detection)
* [Miscellaneous Hi-C tools](#miscellaneous)
* [Papers](#papers)
  * [Methodological Reviews](#methodological-reviews)
  * [General Reviews](#general-reviews)
  * [Normalization](#normalization)
  * [TAD detection](#tad-detection)
  * [TAD prediction](#tad-prediction)
  * [Spectral clustering](#spectral-clustering)


## Pipelines

- `distiller-nf` - Java modular Hi-C mapping pipeline for reproducible data analysis, nextflow pipeline. Alignment, filtering, aggregating Hi-C matrices. https://github.com/mirnylab/distiller-nf

- `Juicer` - Java full pipeline to convert raw reads into Hi-C maps, visualized in Juicebox. Call domains, loops, CTCF binding sites. `.hic` file format for storing multi-resolution Hi-C data. https://github.com/theaidenlab/juicebox/wiki/Download
    - Durand, Neva C., Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. “Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments.” Cell Systems 3, no. 1 (July 2016): 95–98. https://doi.org/10.1016/j.cels.2016.07.002.
    - Rao, Suhas S. P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. “A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping.” Cell 159, no. 7 (December 18, 2014): 1665–80. https://doi.org/10.1016/j.cell.2014.11.021. - Juicer analysis example. TADs defined by frequent interactions. Enriched in CTCF and cohesin members. Five domain types. A1 and A2 enriched in genes. Chr 19 contains 6th pattern B6. Enrichment in different histone modification marks. TADs are preserved across cell types. Yet, differences between Gm12878 and IMR90 were detected. Boundaries detection by scanning image. Refs to the original paper.


- `HiCExplorer` - set of Python scripts to process, normalize, analyze and visualize Hi-C data, Python. https://hicexplorer.readthedocs.io/en/latest/, https://github.com/deeptools/HiCExplorer/

- `hiclib` - Python tools to qc, map, normalize, filter and analyze Hi-C data, https://bitbucket.org/mirnylab/hiclib

- `HiC-bench` - complete pipeline for Hi-C data analysis. https://github.com/NYU-BFX/hic-bench
    - Lazaris, Charalampos, Stephen Kelly, Panagiotis Ntziachristos, Iannis Aifantis, and Aristotelis Tsirigos. “HiC-Bench: Comprehensive and Reproducible Hi-C Data Analysis Designed for Parameter Exploration and Benchmarking.” BMC Genomics 18, no. 1 (December 2017). https://doi.org/10.1186/s12864-016-3387-6.

- `HiC-Pro` - Python and command line-based optimized and flexible pipeline for Hi-C data processing, https://github.com/nservant/HiC-Pro
    - Servant, Nicolas, Nelle Varoquaux, Bryan R. Lajoie, Eric Viara, Chong-Jian Chen, Jean-Philippe Vert, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiC-Pro: An Optimized and Flexible Pipeline for Hi-C Data Processing.” Genome Biology 16 (December 1, 2015): 259. https://doi.org/10.1186/s13059-015-0831-x. - HiC pipeline, references to other pipelines, comparison. From raw reads to normalized matrices. Normalization methods, fast and memory-efficient implementation of iterative correction normalization (ICE). Data format. Using genotyping information to phase contact maps.

- `HiC_Pipeline` - Python-based pipeline performing mapping, filtering, binning, and ICE-correcting Hi-C data, from raw reads (.sra, .fastq) to contact matrices. Additionally, converting to sparse format, performing QC. https://github.com/XiaoTaoWang/HiC_pipeline

- `HiCUP` - Perl-based pipeline, alignment only, output - BAM files. http://www.bioinformatics.babraham.ac.uk/projects/hicup/
    - Wingett, Steven, Philip Ewels, Mayra Furlan-Magaril, Takashi Nagano, Stefan Schoenfelder, Peter Fraser, and Simon Andrews. “HiCUP: Pipeline for Mapping and Processing Hi-C Data.” F1000Research 4 (2015): 1310. https://doi.org/10.12688/f1000research.7334.1. - HiCUP pipeline, alignment only, removes artifacts (religations, duplicate reads) creating BAM files. Details about Hi-C sequencing artefacts. Used in conjunction with other pipelines.

- `HiFi` - Python/C++ tool for extracting restriction fragment resolution Hi-C data. https://github.com/BlanchetteLab/HIFI

- `HiTC` - R package for High Throughput Chromosome Conformation Capture analysis, https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. https://bioconductor.org/packages/release/bioc/html/HiTC.html and https://www.bioconductor.org/packages/devel/bioc/vignettes/HiTC/inst/doc/HiTC.pdf 

- `my5C`- web-based tools, well-documented analysis and visualization of 5S data, http://my5c.umassmed.edu/

- `TADbit` - Python-based pipeline, from iterative mapping, filtering, normalization. Similarity metrics: distance-centric Spearman, first principal eigenvector. TAD detection. TAD boundaries alignment, within 20kb. 3D modeling. Supplementary material - key functions, TAD detection algorithm, boundary comparison. https://github.com/3DGenomes/tadbit
    - Serra, François, Davide Baù, Mike Goodstadt, David Castillo, Guillaume J. Filion, and Marc A. Marti-Renom. “Automatic Analysis and 3D-Modelling of Hi-C Data Using TADbit Reveals Structural Features of the Fly Chromatin Colors.” PLoS Computational Biology 13, no. 7 (July 2017): e1005665. https://doi.org/10.1371/journal.pcbi.1005665.



## Normalization

- `HiCNorm` - removing biases in Hi-C data via Poisson regression, http://www.people.fas.harvard.edu/~junliu/HiCNorm/
    - Hu, Ming, Ke Deng, Siddarth Selvaraj, Zhaohui Qin, Bing Ren, and Jun S. Liu. “HiCNorm: Removing Biases in Hi-C Data via Poisson Regression.” Bioinformatics (Oxford, England) 28, no. 23 (December 1, 2012): 3131–33. https://doi.org/10.1093/bioinformatics/bts570. - Poisson normalization. Also tested negative binomial.

- `HiFive` - handling and normalization or pre-aligned Hi-C and 5C data, https://www.taylorlab.org/software/hifive/
    - Sauria, Michael EG, Jennifer E. Phillips-Cremins, Victor G. Corces, and James Taylor. “HiFive: A Tool Suite for Easy and Efficient HiC and 5C Data Analysis.” Genome Biology 16, no. 1 (December 2015). https://doi.org/10.1186/s13059-015-0806-y. - HiFive - post-processing of aligned Hi-C and 5C data, three normalization approaches: "Binning" - model-based Yaffe & Tanay's method, "Express" - matrix-balancing approach, "Probability" - multiplicative probability model. Judging normalization quality by correlation between matrices. 

- `HiTC` - The HiTC R package was developed to explore high-throughput 'C' data such as 5C or Hi-C. Dedicated R classes as well as standard methods for quality controls, normalization, visualization, and further analysis are also provided. https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization using loess regression on genomic distance, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. 

- `caICB` - chromosome-adjusted iterative correction of biases method. Review of Hi-C data biases and methods to remove them. CNV bias cannot be removed after within-chromosome iterative correction bias methods. The distance-dependent decay of interaction frequencies is modeled by splines. The caICB method aims to minimize the differences across count-distance curves of different chromosomes, chr1 as the reference. A priori knowledge of CNVs is not required. The between- and within-chromosome bias is removed, minimizing the number of significant contacts due to CNVs. https://bitbucket.org/mthjwu/hicapp
    - Wu, Hua-Jun, and Franziska Michor. “A Computational Strategy to Adjust for Copy Number in Tumor Hi-C Data.” Bioinformatics (Oxford, England) 32, no. 24 (15 2016): 3695–3701. https://doi.org/10.1093/bioinformatics/btw540.

- `ICE` - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. 
    - Imakaev, Maxim, Geoffrey Fudenberg, Rachel Patton McCord, Natalia Naumova, Anton Goloborodko, Bryan R. Lajoie, Job Dekker, and Leonid A. Mirny. “Iterative Correction of Hi-C Data Reveals Hallmarks of Chromosome Organization.” Nature Methods 9, no. 10 (October 2012): 999–1003. https://doi.org/10.1038/nmeth.2148. - ICE - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. Assumption - all loci should have equal visibility. Deconvolution into eigenvectors/values. hiclib https://bitbucket.org/mirnylab/hiclib. Good description of the algorithm by Lior Pachter https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/

- `OneD` - CNV bias-correction method, addresses the problem of partial aneuploidy. Bin-centric counts are modeled using negative binomial distribution, and its parameters are estimated using splines. A hidden Markov model is fit to infer copy number for each bin. Each Hi-C matrix entry is corrected by dividing its value by square root of the product of CNVs for the corresponding bins. Reproducibility score (eigenvector decomposition and comparison) to measure improvement in the similarity between replicated Hi-C data. https://github.com/qenvio/dryhic
    - Vidal, Enrique, François le Dily, Javier Quilez, Ralph Stadhouders, Yasmina Cuartero, Thomas Graf, Marc A Marti-Renom, Miguel Beato, and Guillaume J Filion. “OneD: Increasing Reproducibility of Hi-C Samples with Abnormal Karyotypes.” Nucleic Acids Research, January 31, 2018. https://doi.org/10.1093/nar/gky064.


## Reproducibility

- `3DChromatin_ReplicateQC` - Comparison of four Hi-C reproducibility assessment tools, `HiCRep`, `GenomeDISCO`, `HiC-Spector`, `QuASAR-Rep`. Tested the effects of noise, sparsity, resolution. Spearman doesn't work well. All tools performed similarly, worsening expectedly. QuASAR has QC tool measuring the level of noise. https://github.com/kundajelab/3DChromatin_ReplicateQC
    - Yardimci, Galip, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data,” September 14, 2017. doi:10.1101/188755. 

- `HiC-Spector` - reproducibility metric to quantify the similarity between contact maps using spectral decomposition. Decomposing Laplacian matrices and sum the Euclidean distance between eigenvectors. https://github.com/gersteinlab/HiC-spector
    - Yan, Koon-Kiu, Galip Gürkan Yardimci, Chengfei Yan, William S. Noble, and Mark Gerstein. “HiC-Spector: A Matrix Library for Spectral and Reproducibility Analysis of Hi-C Contact Maps.” Bioinformatics (Oxford, England) 33, no. 14 (July 15, 2017): 2199–2201. https://doi.org/10.1093/bioinformatics/btx152.


## Interaction callers

- `Fit-Hi-C` - Python tool for detection of significant chromatin interactions, https://noble.gs.washington.edu/proj/fit-hi-c/
    - Ay, Ferhat, Timothy L. Bailey, and William Stafford Noble. “Statistical Confidence Estimation for Hi-C Data Reveals Regulatory Chromatin Contacts.” Genome Research 24, no. 6 (June 2014): 999–1011. https://doi.org/10.1101/gr.160374.113. - Fit-Hi-C method, Splines to model distance dependence. Model mid-range interaction frequencies, decay with distance. Biases, methods for normalization. Two-step splines - use all dots for first fit, identify and remove outliers, second fit without outliers. Markers of boundaries - insulators, heterochromatin, pluripotent factors. CNVs are enriched in chromatin boundaries. Replication timing data how-to http://www.replicationdomain.com/. Validation Hi-C data. http://chromosome.sdsc.edu/mouse/hi-c/download.html

- `GoTHIC` - R package for peak calling in individual HiC datasets, while accounting for noise. https://www.bioconductor.org/packages/release/bioc/html/GOTHiC.html
    - Mifsud, Borbala, Inigo Martincorena, Elodie Darbo, Robert Sugar, Stefan Schoenfelder, Peter Fraser, and Nicholas M. Luscombe. “GOTHiC, a Probabilistic Model to Resolve Complex Biases and to Identify Real Interactions in Hi-C Data.” Edited by Mark Isalan. PLOS ONE 12, no. 4 (April 5, 2017): e0174744. https://doi.org/10.1371/journal.pone.0174744. - The GOTHiC (genome organisation through HiC) algorithm uses a simple binomial distribution model to simultaneously remove coveralge-associated biases in Hi-C data and detect significant interactions by assuming that the global background interaction frequency of two loci. Use of the Benjamini–Hochberg multiple-testing correction to control for the false discovery rate. 

- `HiCPeaks` - Python CPU-based implementation for BH-FDR and HICCUPS, two peak calling algorithms for Hi-C data, proposed by Rao et al 2014. Text-to-cooler Hi-C data converter, two scripts to call peaks, and one for visualization (creation of a .png file)

- `HOMER` - Perl scripts for normalization, visualization, significant interaction detection, motif discovery. Does not correct for bias. http://homer.ucsd.edu/homer/interactions/


## TAD callers

- `3D-NetMod` - hierarchical, nested, partially overlapping TAD detection using graph theory. Community detection method based on the maximization of network modularity, Louvain-like locally greedy algorithm, repeated several (20) times to avoid local maxima, then getting consensus. Tuning parameters are estimated over sequence search. Benchmarked against TADtree, directionality index, Arrowhead. ICE-normalized data brain data from Geschwind (human data) and Jiang (mouse data) studies. Computationally intensive. Python implementation https://bitbucket.org/creminslab/3dnetmod_method_v1.0_10_06_17
    - Norton, Heidi K., Daniel J. Emerson, Harvey Huang, Jesi Kim, Katelyn R. Titus, Shi Gu, Danielle S. Bassett, and Jennifer E. Phillips-Cremins. “Detecting Hierarchical Genome Folding with Network Modularity.” Nature Methods 15, no. 2 (February 2018): 119–22. https://doi.org/10.1038/nmeth.4560.

- `Armatus` - TAD detection at different resolutions, https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus
    - Filippova, Darya, Rob Patro, Geet Duggal, and Carl Kingsford. “Identification of Alternative Topological Domains in Chromatin.” Algorithms for Molecular Biology 9, no. 1 (2014): 14. doi:10.1186/1748-7188-9-14.https://almob.biomedcentral.com/track/pdf/10.1186/1748-7188-9-14 - Dynamic programming method named “Armatus”. Methods - statistics, intuitive. Consider different resolutions, stability of the detected TADs, persistency across resolution. Only one parameter controlling the size of domains. Their TADs are different from Dixon's TADs. https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus

- `CaTCH` - identification of hierarchical TAD structure, https://github.com/zhanyinx/CaTCH_R
    - Zhan, Yinxiu, Luca Mariani, Iros Barozzi, Edda G. Schulz, Nils Blüthgen, Michael Stadler, Guido Tiana, and Luca Giorgetti. “Reciprocal Insulation Analysis of Hi-C Data Shows That TADs Represent a Functionally but Not Structurally Privileged Scale in the Hierarchical Folding of Chromosomes.” Genome Research 27, no. 3 (2017): 479–90. https://doi.org/10.1101/gr.212803.116. - CaTCH - identification of hierarchical TAD structure. Reciprocal insulation (RI) index. Benchmarked against Dixon's TADs (diTADs). CTCF enrichment as a benchmark, enrichment of TADs in differentially expressed genes. https://github.com/zhanyinx/CaTCH_R

- `cLoops` - loops calling for ChIA-PET, Hi-C and HiChIP. https://github.com/YaqiangCao/cLoops

- `ClusterTAD` - A clustering method for identifying topologically associated domains (TADs) from Hi-C data, https://github.com/BDM-Lab/ClusterTAD
    - Oluwadare, Oluwatosin, and Jianlin Cheng. “ClusterTAD: An Unsupervised Machine Learning Approach to Detecting Topologically Associated Domains of Chromosomes from Hi-C Data.” BMC Bioinformatics 18, no. 1 (November 14, 2017): 480. https://doi.org/10.1186/s12859-017-1931-2. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1931-2 - ClusterTAD paper. Clustering to define TADs. Datasets: simulated Hi-C data with pre-defined TADs https://link.springer.com/article/10.1007%2Fs40484-015-0047-9, RenLab Hi-C and CTCF data http://chromosome.sdsc.edu/mouse/download.html. 

- `EAST` - Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps, https://github.com/ucrbioinfo/EAST
    - Abbas Roayaei Ardakany, Stefano Lonardi, and Marc Herbstritt, “Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps” (Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik GmbH, Wadern/Saarbruecken, Germany, 2017), https://doi.org/10.4230/LIPIcs.WABI.2017.22. - EAST: Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps. Haar-like features (rectangles on images) and a function that quantifies TAD properties: frequency within is high, outside - low, boundaries must be strong. Objective - finding a set of contigious non-overlapping domains maximizing the function. Restricted by maximum length of TADs Text boundaries for enrichment in CTCF, RNP PolII, H3K4me3, H3K27ac. https://github.com/ucrbioinfo/EAST

- `hickit` - TAD calling, phase imputation, 3D modeling and more for diploid single-cell Hi-C (Dip-C) and bulk Hi-C, https://github.com/lh3/hickit

- `HiTAD` - hierarchical TAD identification, different resolutions, correlation with chromosomal compartments, replication timing, gene expression. Adaptive directionality index approach. Data sources, methods for comparing TAD boundaries, reproducibility. H3K4me3 enriched and H3K4me1 depleted at boundaries. TAD boundaries (but not sub-TADs) separate replication timing, A/B compartments, gene expression. https://github.com/XiaoTaoWang/TADLib,  https://pypi.python.org/pypi/TADLib
    - Wang, Xiao-Tao, Wang Cui, and Cheng Peng. “HiTAD: Detecting the Structural and Functional Hierarchies of Topologically Associating Domains from Chromatin Interactions.” Nucleic Acids Research 45, no. 19 (November 2, 2017): e163. https://doi.org/10.1093/nar/gkx735.

- `IC-Finder` - Segmentations of HiC maps into hierarchical interaction compartments, http://membres-timc.imag.fr/Daniel.Jost/DJ-TIMC/Software.html
    - Noelle Haddad, Cedric Vaillant, Daniel Jost. "IC-Finder: inferring robustly the hierarchical organization of chromatin folding." Nucleic Acids Res. 2017 Jun 2; 45(10). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5449546/. TAD-identification procedure that is essentially identical to my idea, including inferring hierarchical structure of TADs based on hierarchical clustering. Article does not explicity discuss A/B compartments so it may still be novel for compartment detection. - TAD finding by combined hierarchical clustering. Correlation-based distance, weighted-mean linkage similarity. Variance criterion to define cluster boundaries. http://membres-timc.imag.fr/Daniel.Jost/DJ-TIMC/Software.html

- Insulation score: Giorgetti, Luca, Bryan R. Lajoie, Ava C. Carter, Mikael Attia, Ye Zhan, Jin Xu, Chong Jian Chen, et al. “Structural Organization of the Inactive X Chromosome in the Mouse.” Nature 535, no. 7613 (28 2016): 575–79. https://doi.org/10.1038/nature18589. - https://github.com/dekkerlab/cworld-dekker/tree/master/scripts/perl matrix2insulation.pl, Parameters: -is 480000 -ids 320000 -im iqrMean -nt 0 -ss 160000 -yb 1.5 -nt 0 -bmoe 0. 
    - Used in Yardımcı, Galip Gürkan, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data.” BioRxiv, January 1, 2018. https://doi.org/10.1101/188755. - TADs detected by insulation score are robust to resolution and noise

- `HiCseg` - TAD detection by maximization of likelihood based block-wise segmentation model, https://cran.r-project.org/web/packages/HiCseg/index.html
    - Lévy-Leduc, Celine, M. Delattre, T. Mary-Huard, and S. Robin. “Two-Dimensional Segmentation for Analyzing Hi-C Data.” Bioinformatics (Oxford, England) 30, no. 17 (September 1, 2014): i386-392. https://doi.org/10.1093/bioinformatics/btu443. - HiCseg paper. TAD detection by maximization of likelihood based block-wise segmentation model, HiCseg R package. 2D segmentation rephrased as 1D segmentation - not contours, but borders. Statistical framework, solved with dynamic programming. Dixon data as gold standard. Hausdorff distance to compare segmentation quality, https://en.wikipedia.org/wiki/Hausdorff_distance. Parameters (from TopDom paper): nb_change_max = 500, distrib = 'G' and model = 'Dplus'.

- `OnTAD` - hierarchical TAD caller, Optimal Nested TAD caller. Sliding window, adaptive local minimum search algorithm, similar to TOPDOM. https://github.com/anlin00007/OnTAD
    - An, Lin, Tao Yang, Jiahao Yang, Johannes Nuebler, Qunhua Li, and Yu Zhang. “Hierarchical Domain Structure Reveals the Divergence of Activity among TADs and Boundaries,” July 3, 2018. https://doi.org/10.1101/361147. - Intro about TADs, Dixon's directionality index, Insulation score. Other hierarchical callers - TADtree, rGMAP, Arrowhead, 3D-Net, IC-Finder. Limitations of current callers - ad hoc thresholds, sensitivity to sequencing depth and mapping resolution, long running time and large memory usage, insufficient performance evaluation. Boundaries are asymmetric - some have more contacts with other boundaries, support for asymmetric loop extrusion model. Performance comparison with DomainCaller, rGMAP, Arrowhead, TADtree. Stronger enrichment of CTCF and two cohesin proteins RAD21 and SMC3. TAD-adjR^2 metric quantifying proportion of variance in the contact frequencies explained by TAD moundaries. Reproducibility of TAD boundaries - Jaccard index, tested at different sequencing depths and resolutions. Boundaries of hierarchical TADs are more active - more CTCF, epigenomic features, TFBSs expressed genes. Super-boundaries - shared by 5 or more TADs, highly active. Rao-Huntley 2014 Gm12878 data. Distance correction - subtracting the mean counts at each distance.

- `TADtree` - TADtree is an algorithm the identification of hierarchical topological domains in Hi-C data, http://compbio.cs.brown.edu/software/
    - Weinreb, Caleb, and Benjamin J. Raphael. “Identification of Hierarchical Chromatin Domains.” Bioinformatics (Oxford, England) 32, no. 11 (June 1, 2016): 1601–9. https://doi.org/10.1093/bioinformatics/btv485. - TADtree paper. Hierarchical (nested) TAD identification. Two ways of TAD definition: 1D and 2D. Normalization by distance. Enrichment over background. Deep statistics of the method. How to compare TADs (VI measure (`vi.dist` in https://cran.r-project.org/web/packages/mcclust/mcclust.pdf), Precision/recall using Dixon as the true set, Fig. 5: number of TADs, TAD size boxplots, Enrichment within 50kb of a TAD boundary - CTCF, PolII, H3K4me3, housekeeping genes - stronger enrichment the better). http://compbio.cs.brown.edu/software/

- `TopDom` - An efficient and Deterministic Method for identifying Topological Domains in Genomes, http://zhoulab.usc.edu/TopDom/
    - Shin, Hanjun, Yi Shi, Chao Dai, Harianto Tjong, Ke Gong, Frank Alber, and Xianghong Jasmine Zhou. “TopDom: An Efficient and Deterministic Method for Identifying Topological Domains in Genomes.” Nucleic Acids Research 44, no. 7 (April 20, 2016): e70. https://doi.org/10.1093/nar/gkv1505. - TopDom paper. Review of other methods. Method is based on general observation that within-TAD interactions are stronger than between-TAD. binSignal value as the average of nearby contact frequency, fitting a curve, finding local minima, test them for significance. Fast, takes linear time. Detects similar domains to HiCseq and Dixon's directionaliry index. Found expected enrichment in CTCF, histone marks. Housekeeping genes and overall gene density are close to TAD boundaries, differentially expressed genes are not. Figure 7 - how to detect common/unique boundaries using Jaccard-like statistics. http://zhoulab.usc.edu/TopDom/


## Prediction of 3D features

- `3DEpiLoop` - prediction of 3D interactions from 1D epigenomic profiles using Random Forest trained on CTCF peaks (histone modifications are the most important predictors, and TFBSs). https://bitbucket.org/4dnucleome/3depiloop
    - Al Bkhetan, Ziad, and Dariusz Plewczynski. “Three-Dimensional Epigenome Statistical Model: Genome-Wide Chromatin Looping Prediction.” Scientific Reports 8, no. 1 (December 2018). https://doi.org/10.1038/s41598-018-23276-8.


## SNP-oriented

- `iRegNet3D` - Integrated Regulatory Network 3D (iRegNet3D) is a high-resolution regulatory network comprised of interfaces of all known transcription factor (TF)-TF, TF-DNA interaction interfaces, as well as chromatin-chromatin interactions and topologically associating domain (TAD) information from different cell lines.  
    - Goal: SNP interpretation
    - Input: One or several SNPs, rsIDs or genomic coordinates.
    - Output: For one or two SNPs, on-screen information of their disease-related info, connection over TF-TF and chromatin interaction networks, and whether they interact in 3D and located within TADs. For multiple SNPs, same info downloadable as text files.

- `3DSNP` - A database linking noncoding SNPs to 3D interacting genes. http://cbportal.org/3dsnp/
    - Lu, Yiming, Cheng Quan, Hebing Chen, Xiaochen Bo, and Chenggang Zhang. “3DSNP: A Database for Linking Human Noncoding SNPs to Their Three-Dimensional Interacting Genes.” Nucleic Acids Research 45, no. D1 (January 4, 2017): D643–49. https://doi.org/10.1093/nar/gkw1022. - 3DSNP database integrating SNP epigenomic annotations with chromatin loops. Linear closest gene, 3D interacting gene, eQTL, 3D interacting SNP, chromatin states, TFBSs, conservation. For individual SNPs.

- HUGIn, tissue-specific Hi-C linear display of anchor position and around. Overlay gene expression and epigenomic data. Association of SNPs with genes based on Hi-C interactions. Tissue-specific. http://yunliweb.its.unc.edu/HUGIn/
    - Martin, Joshua S, Zheng Xu, Alex P Reiner, Karen L Mohlke, Patrick Sullivan, Bing Ren, Ming Hu, and Yun Li. “HUGIn: Hi-C Unifying Genomic Interrogator.” BioRxiv, 2017, 117531.


## Structural variant detection

- `hic_breakfinder` - SV identification in Hi-C data. https://github.com/dixonlab/hic_breakfinder
    - Dixon, Jesse R., Jie Xu, Vishnu Dileep, Ye Zhan, Fan Song, Victoria T. Le, Galip Gürkan Yardımcı, et al. “Integrative Detection and Analysis of Structural Variation in Cancer Genomes.” Nature Genetics, September 10, 2018. https://doi.org/10.1038/s41588-018-0195-8. - Detection of structural variants (SV) by integrating optical mapping, Hi-C, and WGS. Custom pipeline using LUMPY, Delly, Control-FREEC software. New Hi-C data on 14 cancer cell lines and 21 previously published datasets. Integration of the detected SVs with genomic annotations, including replication timing. Supplementary data with SVs resolved by individual methods and integrative approaches.


## Miscellaneous

- `HiCPlus` - increasing resolution of Hi-C data using convolutional neural network. Basically, smoothing parts of Hi-C image, then binning into smaller parts. Performs better than bilinear/biqubic smoothing. https://github.com/zhangyan32/HiCPlus
    - Zhang, Yan, Lin An, Ming Hu, Jijun Tang, and Feng Yue. “HiCPlus: Resolution Enhancement of Hi-C Interaction Heatmap,” March 1, 2017. https://doi.org/10.1101/112631.




## Papers

### Methodological Reviews

- Nicoletti, Chiara, Mattia Forcato, and Silvio Bicciato. “Computational Methods for Analyzing Genome-Wide Chromosome Conformation Capture Data.” Current Opinion in Biotechnology 54 (December 2018): 98–105. https://doi.org/10.1016/j.copbio.2018.01.023. - 3C-Hi-C tools review, Table 1 lists categorizes main tools, Figure 1 displays all steps in technology and analysis (alignment, resolution, normalization, including accounting for CNVs, A/B compartments, TAD detection, visualization). Concise description of all tools.

- Waldispühl, Jérôme, Eric Zhang, Alexander Butyaev, Elena Nazarova, and Yan Cyr. “Storage, Visualization, and Navigation of 3D Genomics Data.” Methods, May 2018. https://doi.org/10.1016/j.ymeth.2018.05.008. - Review of tools for visualization of 3C-Hi-C data, challenges, analysis (Table 1). Data formats (hic, cool, BUTLR, ccmap). Database to quickly access 3D data. Details of each visualization tool in Section 4

- Forcato, Mattia, Chiara Nicoletti, Koustav Pal, Carmen Maria Livi, Francesco Ferrari, and Silvio Bicciato. “Comparison of Computational Methods for Hi-C Data Analysis.” Nature Methods, June 12, 2017. https://doi.org/10.1038/nmeth.4325. - Hi-C processing and TAD calling tools benchmarking, [Table 1](https://www.nature.com/articles/nmeth.4325/tables/1), simulated (Lun and Smyth method) and real data. Notes about pluses and minuses of each tool. TAD reproducibility is higher than chromatin interactions, increases with larger number of reads. Consistent enrichment of TAD boundaries in CTCF, irrespectively of TAD caller. Hi-C replication is poor, just a bit more than random. Supplementary table 2 - technical details about each program, Supplementary Note 1 - Hi-C preprocessing tools, Supplementary Note 2 - TAD callers. Supplementary note 3 - how to simulate Hi-C data. Supplementary note 6 - how to install tools. https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf

- Yardımcı, Galip Gürkan, and William Stafford Noble. “Software Tools for Visualizing Hi-C Data.” Genome Biology 18, no. 1 (December 2017). https://doi.org/10.1186/s13059-017-1161-y. - Hi-C technology, data, and visualization review. Suggestion about graph representation.

- Ay, Ferhat, and William S. Noble. “Analysis Methods for Studying the 3D Architecture of the Genome.” Genome Biology 16 (September 2, 2015): 183. https://doi.org/10.1186/s13059-015-0745-7. - Hi-C technology and methods review. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7#Tab1) - list of tools. Biases, normalization, matrix balancing. Extracting significant contacts, obs/exp ratio, parametric (powerlaw, neg binomial, double exponential), non-parametric (splines). 3D enrichment. References. TAD identification, directionality index. Outlook, importance of comparative analysis

### General Reviews

- Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation.” Molecular Systems Biology 11, no. 12 (December 23, 2015): 852–852. doi:10.15252/msb.20156492. http://msb.embopress.org/content/msb/11/12/852.full.pdf - 3D genome organization parts. Well-written and detailed. References. Technologies: FISH, 3C. 4C, 5C, Hi-C, GCC, TCC, ChIA-PET. Typical resolution - 40bp to 1Mb. LADs - conserved, but some are cell type-specific. Chromosome territories. Cell type-specific. inter-chromosomal interactions may be important to define cell-specific interactions. A/B compartments identified by PCA. Chromatin loops, marked by CTCF and Cohesin binding, sometimes, with Mediator. Transcription factories

- Dekker, Job, Marc A. Marti-Renom, and Leonid A. Mirny. “Exploring the Three-Dimensional Organization of Genomes: Interpreting Chromatin Interaction Data.” Nature Reviews. Genetics 14, no. 6 (June 2013): 390–403. https://doi.org/10.1038/nrg3454. https://www.nature.com/articles/nrg3454 - 3D genome review. Chromosomal territories, transcription factories. Details of each 3C technology. Exponential decay of interaction frequencies. Box 2: A/B compartments (several Mb), TAD definition, size (hundreds of kb). TADs are largely stable, A/B compartments are tissue-specific. Adjacent TADs are not necessarily of opposing signs, may jointly form A/B compartments. Genes co-expression, enhancer-promoters interactions are confined to TADs. 3D modeling.

- Witten, Daniela M., and William Stafford Noble. “On the Assessment of Statistical Significance of Three-Dimensional Colocalization of Sets of Genomic Elements.” Nucleic Acids Research 40, no. 9 (May 2012): 3849–55. https://doi.org/10.1093/nar/gks012.

### Normalization

- Yaffe, Eitan, and Amos Tanay. “Probabilistic Modeling of Hi-C Contact Maps Eliminates Systematic Biases to Characterize Global Chromosomal Architecture.” Nature Genetics 43, no. 11 (November 2011): 1059–65. https://doi.org/10.1038/ng.947. - Sources of biases: 1) non-specific ligation (large distance between pairs); 2) length of each ligated fragments; 3) CG content and nucleotide composition; 4) Mappability. Normalization. Enrichment of long-range interactions in active promoters. General aggregation of active chromosomal domains. Chromosomal territories, high-activity and two low-activity genomic clusters

### TAD detection

- Dali, Rola, and Mathieu Blanchette. “A Critical Assessment of Topologically Associating Domain Prediction Tools.” Nucleic Acids Research 45, no. 6 (April 7, 2017): 2994–3005. doi:10.1093/nar/gkx145. - TAD definition, tools. Meta-TADs, hierarchy, overlapping TADs. HiCPlotter for visualization. Manual annotation as a gold standard. Sequencing depth and resolution affects things. Code, manual annotations

- Olivares-Chauvet, Pedro, Zohar Mukamel, Aviezer Lifshitz, Omer Schwartzman, Noa Oded Elkayam, Yaniv Lubling, Gintaras Deikus, Robert P. Sebra, and Amos Tanay. “Capturing Pairwise and Multi-Way Chromosomal Conformations Using Chromosomal Walks.” Nature 540, no. 7632 (November 30, 2016): 296–300. https://doi.org/10.1038/nature20158. - TADs organize chromosomal territories. Active and inactive TAD properties. Methods: Good mathematical description of insulation score calculations. Filter TADs smaller than 250kb. Inter-chromosomal contacts are rare, ~7-10%. Concatemers (more than two contacts) are unlikely.


- Crane, Emily, Qian Bian, Rachel Patton McCord, Bryan R. Lajoie, Bayly S. Wheeler, Edward J. Ralston, Satoru Uzawa, Job Dekker, and Barbara J. Meyer. “Condensin-Driven Remodelling of X Chromosome Topology during Dosage Compensation.” Nature 523, no. 7559 (July 9, 2015): 240–44. https://doi.org/10.1038/nature14450. - InsulationScore, https://github.com/dekkerlab/crane-nature-2015 - Insulation score to define TADs - sliding square along the diagonal, aggregating signal within it. This aggregated score is normalized, and binned into TADs, boundaries. See Methods and implementation at https://github.com/dekkerlab/crane-nature-2015. ICE normalized data. OK to analyze data at two different resolutions

- "Hierarchical Regulatory Domain Inference from Hi-C Data" - presentation by Bartek Wilczyński about TAD detection, existing algorithms, new SHERPA and OPPA methods. [Video](https://simons.berkeley.edu/talks/bartek-wilczynski-03-10-16), [PDF](https://simons.berkeley.edu/sites/default/files/docs/4588/2016-03-10-simons-institute-wilczynski.pdf), [Web site](http://regulomics.mimuw.edu.pl/wp/), [GitHub](https://github.com/regulomics/) - SHERPA and OPPA code there.

### TAD prediction

- Bednarz, Paweł, and Bartek Wilczyński. “Supervised Learning Method for Predicting Chromatin Boundary Associated Insulator Elements.” Journal of Bioinformatics and Computational Biology 12, no. 06 (December 2014): 1442006. doi:10.1142/S0219720014420062. http://www.worldscientific.com/doi/pdf/10.1142/S0219720014420062 - Predicting TAD boundaries using training data, and making new predictions. Bayesian network (BNFinder method), random forest vs. basic k-means clustering, ChromHMM, cdBEST. Using sequence k-mers and ChIP-seq data from modENCODE for prediction - CTCF ChIP-seq performs best. Used Boruta package for feature selection. Bayesian network performs best. To read on their BNFinder method

### Spectral clustering

- Y. X Rachel Wang, Purnamrita Sarkar, Oana Ursu, Anshul Kundaje and Peter J. Bickel, "Network modelling of topological domains using Hi-C data", https://arxiv.org/abs/1707.09587.  - TAD analysis using graph theoretical (network-based) methods. Treats TADs as a "community" within the network. Shows that naive spectral clustering is generally ineffective, leaving gaps in the data. 

- Liu, Sijia, Pin-Yu Chen, Alfred Hero, and Indika Rajapakse. “Dynamic Network Analysis of the 4D Nucleome.” BioRxiv, January 1, 2018. https://doi.org/10.1101/268318. - Temporal Hi-C data analysis using graph theory. Integrated with RNA-seq data. Network-based approaches such as von Neumann graph entropy, network centrality, and multilayer network theory are applied to reveal universal patterns of the dynamic genome. Toeplitz normalization. Graph Laplasian matrix. Detailed statistics.

- Norton, Heidi K, Harvey Huang, Daniel J Emerson, Jesi Kim, Shi Gu, Danielle S Bassett, and Jennifer E Phillips-Cremins. “Detecting Hierarchical 3-D Genome Domain Reconfiguration with Network Modularity,” November 22, 2016. https://doi.org/10.1101/089011. - Graph theory for TAD identification. Louvain-like local greedy algorithm to maximize network modularity. Vary resolution parameter, hierarchical TAD identification. Hierarchical spatial variance minimization method. ROC analysis to quantify performance. Adjusted RAND score to quantify TAD overlap.

- Chen, Jie, Alfred O. Hero, and Indika Rajapakse. “Spectral Identification of Topological Domains.” Bioinformatics (Oxford, England) 32, no. 14 (15 2016): 2151–58. https://doi.org/10.1093/bioinformatics/btw221. - Spectral algorithm to define TADs. Laplacian graph segmentation using Fiedler vector iteratively. Toeplitz normalization to remove distance effect. Spectral TADs do not overlap with Dixon's, but better overlap with CTCF.

- Fotuhi Siahpirani, Alireza, Ferhat Ay, and Sushmita Roy. “A Multi-Task Graph-Clustering Approach for Chromosome Conformation Capture Data Sets Identifies Conserved Modules of Chromosomal Interactions.” Genome Biology 17, no. 1 (December 2016). https://doi.org/10.1186/s13059-016-0962-8 - Arboretum-Hi-C - a multitask spectral clustering method to identify differences in genomic architecture. Intro about the 3D genome organization, TAD differences and conservation. Assessment of different clustering approaches using different distance measures, as well as raw contacts. Judging clustering quality by enrichment in regulatory genomic signals (Histone marks, LADs, early vs. late replication timing, TFs like POLII, TAF, TBP, CTCF, P300, CMYC, cohesin components, LADs, replication timing, SINE, LINE, LTR) and by numerical methods (Davies-Bouldin index, silhouette score, others). Although spectral clustering on contact counts performed best, spectral + Spearman correlation was chosen. Comparing cell types identifies biologically relevant differences as quantified by enrichment. Peak counts or average signal within regions were used for enrichment. Data https://zenodo.org/record/49767, and Arboretum-HiC https://bitbucket.org/roygroup/arboretum-hic



