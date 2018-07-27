# TAD detection

- `Armatus` - TAD detection at different resolutions, https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus
    - Filippova, Darya, Rob Patro, Geet Duggal, and Carl Kingsford. “Identification of Alternative Topological Domains in Chromatin.” Algorithms for Molecular Biology 9, no. 1 (2014): 14. doi:10.1186/1748-7188-9-14.https://almob.biomedcentral.com/track/pdf/10.1186/1748-7188-9-14 - Dynamic programming method named “Armatus”. Methods - statistics, intuitive. Consider different resolutions, stability of the detected TADs, persistency across resolution. Only one parameter controlling the size of domains. Their TADs are different from Dixon's TADs. https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus

- `CaTCH` - identification of hierarchical TAD structure, https://github.com/zhanyinx/CaTCH_R
    - Zhan, Yinxiu, Luca Mariani, Iros Barozzi, Edda G. Schulz, Nils Blüthgen, Michael Stadler, Guido Tiana, and Luca Giorgetti. “Reciprocal Insulation Analysis of Hi-C Data Shows That TADs Represent a Functionally but Not Structurally Privileged Scale in the Hierarchical Folding of Chromosomes.” Genome Research 27, no. 3 (2017): 479–90. https://doi.org/10.1101/gr.212803.116. - CaTCH - identification of hierarchical TAD structure. Reciprocal insulation (RI) index. Benchmarked against Dixon's TADs (diTADs). CTCF enrichment as a benchmark, enrichment of TADs in differentially expressed genes. https://github.com/zhanyinx/CaTCH_R

- `cLoops` - loops calling for ChIA-PET, Hi-C and HiChIP. https://github.com/YaqiangCao/cLoops

- `ClusterTAD` - A clustering method for identifying topologically associated domains (TADs) from Hi-C data, https://github.com/BDM-Lab/ClusterTAD
    - Oluwadare, Oluwatosin, and Jianlin Cheng. “ClusterTAD: An Unsupervised Machine Learning Approach to Detecting Topologically Associated Domains of Chromosomes from Hi-C Data.” BMC Bioinformatics 18, no. 1 (November 14, 2017): 480. https://doi.org/10.1186/s12859-017-1931-2. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1931-2 - ClusterTAD paper. Clustering to define TADs. Datasets: simulated Hi-C data with pre-defined TADs https://link.springer.com/article/10.1007%2Fs40484-015-0047-9, RenLab Hi-C and CTCF data http://chromosome.sdsc.edu/mouse/download.html. 

- `EAST` - Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps, https://github.com/ucrbioinfo/EAST
    - Abbas Roayaei Ardakany, Stefano Lonardi, and Marc Herbstritt, “Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps” (Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik GmbH, Wadern/Saarbruecken, Germany, 2017), https://doi.org/10.4230/LIPIcs.WABI.2017.22. - EAST: Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps. Haar-like features (rectangles on images) and a function that quantifies TAD properties: frequency within is high, outside - low, boundaries must be strong. Objective - finding a set of contigious non-overlapping domains maximizing the function. Restricted by maximum length of TADs Text boundaries for enrichment in CTCF, RNP PolII, H3K4me3, H3K27ac. https://github.com/ucrbioinfo/EAST

- `IC-Finder` - Segmentations of HiC maps into hierarchical interaction compartments, http://membres-timc.imag.fr/Daniel.Jost/DJ-TIMC/Software.html
    - Noelle Haddad, Cedric Vaillant, Daniel Jost. "IC-Finder: inferring robustly the hierarchical organization of chromatin folding." Nucleic Acids Res. 2017 Jun 2; 45(10). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5449546/. TAD-identification procedure that is essentially identical to my idea, including inferring hierarchical structure of TADs based on hierarchical clustering. Article does not explicity discuss A/B compartments so it may still be novel for compartment detection. - TAD finding by combined hierarchical clustering. Correlation-based distance, weighted-mean linkage similarity. Variance criterion to define cluster boundaries. http://membres-timc.imag.fr/Daniel.Jost/DJ-TIMC/Software.html

- Insulation score: Giorgetti, Luca, Bryan R. Lajoie, Ava C. Carter, Mikael Attia, Ye Zhan, Jin Xu, Chong Jian Chen, et al. “Structural Organization of the Inactive X Chromosome in the Mouse.” Nature 535, no. 7613 (28 2016): 575–79. https://doi.org/10.1038/nature18589. - https://github.com/dekkerlab/cworld-dekker/tree/master/scripts/perl matrix2insulation.pl, Parameters: -is 480000 -ids 320000 -im iqrMean -nt 0 -ss 160000 -yb 1.5 -nt 0 -bmoe 0. 
    - Used in Yardımcı, Galip Gürkan, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data.” BioRxiv, January 1, 2018. https://doi.org/10.1101/188755. - TADs detected by insulation score are robust to resolution and noise

- `HiCseg` - TAD detection by maximization of likelihood based block-wise segmentation model, https://cran.r-project.org/web/packages/HiCseg/index.html
    - Lévy-Leduc, Celine, M. Delattre, T. Mary-Huard, and S. Robin. “Two-Dimensional Segmentation for Analyzing Hi-C Data.” Bioinformatics (Oxford, England) 30, no. 17 (September 1, 2014): i386-392. https://doi.org/10.1093/bioinformatics/btu443. - HiCseg paper. TAD detection by maximization of likelihood based block-wise segmentation model, HiCseg R package. 2D segmentation rephrased as 1D segmentation - not contours, but borders. Statistical framework, solved with dynamic programming. Dixon data as gold standard. Hausdorff distance to compare segmentation quality, https://en.wikipedia.org/wiki/Hausdorff_distance. Parameters (from TopDom paper): nb_change_max = 500, distrib = 'G' and model = 'Dplus'.

- `TADtree` - TADtree is an algorithm the identification of hierarchical topological domains in Hi-C data, http://compbio.cs.brown.edu/software/
    - Weinreb, Caleb, and Benjamin J. Raphael. “Identification of Hierarchical Chromatin Domains.” Bioinformatics (Oxford, England) 32, no. 11 (June 1, 2016): 1601–9. https://doi.org/10.1093/bioinformatics/btv485. - TADtree paper. Hierarchical (nested) TAD identification. Two ways of TAD definition: 1D and 2D. Normalization by distance. Enrichment over background. Deep statistics of the method. How to compare TADs (VI measure (`vi.dist` in https://cran.r-project.org/web/packages/mcclust/mcclust.pdf), Precision/recall using Dixon as the true set, Fig. 5: number of TADs, TAD size boxplots, Enrichment within 50kb of a TAD boundary - CTCF, PolII, H3K4me3, housekeeping genes - stronger enrichment the better). http://compbio.cs.brown.edu/software/

- `TopDom` - An efficient and Deterministic Method for identifying Topological Domains in Genomes, http://zhoulab.usc.edu/TopDom/
    - Shin, Hanjun, Yi Shi, Chao Dai, Harianto Tjong, Ke Gong, Frank Alber, and Xianghong Jasmine Zhou. “TopDom: An Efficient and Deterministic Method for Identifying Topological Domains in Genomes.” Nucleic Acids Research 44, no. 7 (April 20, 2016): e70. https://doi.org/10.1093/nar/gkv1505. - TopDom paper. Review of other methods. Method is based on general observation that within-TAD interactions are stronger than between-TAD. binSignal value as the average of nearby contact frequency, fitting a curve, finding local minima, test them for significance. Fast, takes linear time. Detects similar domains to HiCseq and Dixon's directionaliry index. Found expected enrichment in CTCF, histone marks. Housekeeping genes and overall gene density are close to TAD boundaries, differentially expressed genes are not. Figure 7 - how to detect common/unique boundaries using Jaccard-like statistics. http://zhoulab.usc.edu/TopDom/




## Papers

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













