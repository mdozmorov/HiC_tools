# A collection of tools and papers related to Hi-C data analysis

Tools in each section are being resorted newest on top (previously, alphabetically). Related repositories: https://github.com/mdozmorov/HiC_data, https://github.com/mdozmorov/scHiC_notes. Issues and/or Pull requests to add other data are welcome!

# Table of content

* [Pipelines for Hi-C data processing](#pipelines)
  * [Mirnylab tools](#mirnylab-tools)
  * [Capture-C](#capture-c)
  * [HiChIP](#hichip)
  * [4C](#4c)
  * [CUT&RUN](#cut-run)
* [Resolution improvement](#resolution-improvement)
* [Normalization of Hi-C data](#normalization)
  * [CNV-aware normalization](#cnv-aware-normalization)
* [Reproducibility and QC of Hi-C data](#reproducibility)
* [Significant interaction (peak) callers](#significant-interaction-peak-callers)
* [Differential interactions](#differential-interactions)
* [TAD callers](#tad-callers)
* [Prediction of 3D features](#prediction-of-3d-features)
* [SNP-oriented Hi-C analysis](#snp-oriented)
* [CNV and Structural variant detection](#cnv-and-structural-variant-detection)
* [Visualization](#visualization)
* [De novo genome scaffolding](#de-novo-genome-scaffolding)
* [3D reconstruction](#3d-reconstruction)
* [Papers](#papers)
  * [Methodological Reviews](#methodological-reviews)
  * [General Reviews](#general-reviews)
  * [Normalization](#normalization)
  * [TAD detection](#tad-detection)
  * [TAD prediction](#tad-prediction)
  * [TAD dynamics](#tad-dynamics)
  * [Spectral clustering](#spectral-clustering)
* [URLs](#urls)

## Pipelines

- A list of available pipelines, URLs. [pipelines_list.csv](pipelines_list.csv), [Source](https://link.springer.com/protocol/10.1007%2F978-1-4939-8766-5_16)
- Available analysis options in each pipeline. [pipeline_comparison.csv](pipeline_comparison.csv), [Source](https://link.springer.com/protocol/10.1007%2F978-1-4939-8766-5_16)
- [Table summarizing functionality of Hi-C data analysis tools](https://www.sciencedirect.com/science/article/pii/S1672022918304339?via%3Dihub#t0005)
- Review of Hi-C, Capture-C, and Capture-C technologies, their computational preprocessing. Experimental protocols, similarities and differences, types of reads (figures), details of alignment, read orientation, elimination of artefacts, quality metrics. Brief overview of preprocessing tools. Example preprocessing of three types of data. Java tool for preprocessing all types of data, Diachromatic (Differential Analysis of Chromatin Interactions by Capture),  https://github.com/TheJacksonLaboratory/diachromatic, GOPHER (Generator Of Probes for capture Hi-C Experiments at high Resolution) for genome cutting, probe design,  https://github.com/TheJacksonLaboratory/Gopher
    - Hansen, Peter, Michael Gargano, Jochen Hecht, Jonas Ibn-Salem, Guy Karlebach, Johannes T. Roehr, and Peter N. Robinson. “Computational Processing and Quality Control of Hi-C, Capture Hi-C and Capture-C Data.” Genes 10, no. 7 (July 18, 2019): 548. https://doi.org/10.3390/genes10070548.
- Workshop on measuring, analyzing, and visualizing the 3D genome with Hi-C data. Presentations (PDFs, PPTX) and Jupyter notebooks. Cooler, HiCGlass, HiPlier. https://github.com/hms-dbmi/hic-data-analysis-bootcamp

- `cword` - perl cworld module and collection of utility/analysis scripts for C data (3C, 4C, 5C, Hi-C). https://github.com/dekkerlab/cworld-dekker

- `Juicer` - Java full pipeline to convert raw reads into Hi-C maps, visualized in Juicebox. Call domains, loops, CTCF binding sites. `.hic` file format for storing multi-resolution Hi-C data. https://github.com/theaidenlab/juicebox/wiki/Download
    - Durand, Neva C., Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. “Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments.” Cell Systems 3, no. 1 (July 2016): 95–98. https://doi.org/10.1016/j.cels.2016.07.002.
    - Rao, Suhas S. P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. “A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping.” Cell 159, no. 7 (December 18, 2014): 1665–80. https://doi.org/10.1016/j.cell.2014.11.021. - Juicer analysis example. TADs defined by frequent interactions. Enriched in CTCF and cohesin members. Five domain types. A1 and A2 enriched in genes. Chr 19 contains 6th pattern B6. Enrichment in different histone modification marks. TADs are preserved across cell types. Yet, differences between Gm12878 and IMR90 were detected. Boundaries detection by scanning image. Refs to the original paper.

- `GITAR` (HiCtool) - full Hi-C pre-processing, normalization, TAD detection, and visualization. Python scripts wrapping other tools.Table 1 summarizes functionality of existing tools. https://www.genomegitar.org/https://github.com/Zhong-Lab-UCSD/HiCtool
    - Calandrelli, Riccardo, Qiuyang Wu, Jihong Guan, and Sheng Zhong. “GITAR: An Open Source Tool for Analysis and Visualization of Hi-C Data.” Genomics, Proteomics & Bioinformatics 16, no. 5 (2018): 365–72. https://doi.org/10.1016/j.gpb.2018.06.006.

- `HiCExplorer` - set of Python scripts to process, normalize, analyze and visualize Hi-C data, Python. https://hicexplorer.readthedocs.io/en/latest/, https://github.com/deeptools/HiCExplorer/

- `HiCdat` - Hi-C processing pipeline and downstream analysis/visualization. Analyses: normalization, correlation, visualization, comparison, distance decay, PCA, interaction enrichment test, epigenomic enrichment/depletion. https://github.com/MWSchmid/HiCdat
    - Schmid, Marc W., Stefan Grob, and Ueli Grossniklaus. “HiCdat: A Fast and Easy-to-Use Hi-C Data Analysis Tool.” BMC Bioinformatics 16 (September 3, 2015): 277. https://doi.org/10.1186/s12859-015-0678-x.

- `HiCpipe` - an efficient Hi-C data processing pipeline. It is based on Juicer and HiC-pro which combines the advatages of these two processing pipelines. HiCpipe is much faster than Juicer and HiC-pro and can output multile features of Hi-C maps. https://github.com/ChenFengling/HiCpipe

- `HiC-bench` - complete pipeline for Hi-C data analysis. https://github.com/NYU-BFX/hic-bench
    - Lazaris, Charalampos, Stephen Kelly, Panagiotis Ntziachristos, Iannis Aifantis, and Aristotelis Tsirigos. “HiC-Bench: Comprehensive and Reproducible Hi-C Data Analysis Designed for Parameter Exploration and Benchmarking.” BMC Genomics 18, no. 1 (December 2017). https://doi.org/10.1186/s12864-016-3387-6.

- `HiC-Pro` - Python and command line-based optimized and flexible pipeline for Hi-C data processing, https://github.com/nservant/HiC-Pro
    - Servant, Nicolas, Nelle Varoquaux, Bryan R. Lajoie, Eric Viara, Chong-Jian Chen, Jean-Philippe Vert, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiC-Pro: An Optimized and Flexible Pipeline for Hi-C Data Processing.” Genome Biology 16 (December 1, 2015): 259. https://doi.org/10.1186/s13059-015-0831-x. - HiC pipeline, references to other pipelines, comparison. From raw reads to normalized matrices. Normalization methods, fast and memory-efficient implementation of iterative correction normalization (ICE). Data format. Using genotyping information to phase contact maps.

- `HiC_Pipeline` - Python-based pipeline performing mapping, filtering, binning, and ICE-correcting Hi-C data, from raw reads (.sra, .fastq) to contact matrices. Additionally, converting to sparse format, performing QC. https://github.com/XiaoTaoWang/HiC_pipeline

- `HiCUP` - Perl-based pipeline, alignment only, output - BAM files. http://www.bioinformatics.babraham.ac.uk/projects/hicup/
    - Wingett, Steven, Philip Ewels, Mayra Furlan-Magaril, Takashi Nagano, Stefan Schoenfelder, Peter Fraser, and Simon Andrews. “HiCUP: Pipeline for Mapping and Processing Hi-C Data.” F1000Research 4 (2015): 1310. https://doi.org/10.12688/f1000research.7334.1. - HiCUP pipeline, alignment only, removes artifacts (religations, duplicate reads) creating BAM files. Details about Hi-C sequencing artefacts. Used in conjunction with other pipelines.

- `HiTC` - R package for High Throughput Chromosome Conformation Capture analysis, https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. https://bioconductor.org/packages/release/bioc/html/HiTC.html and https://www.bioconductor.org/packages/devel/bioc/vignettes/HiTC/inst/doc/HiTC.pdf 

- `mHi-C` - recovering alignment of multi-mapped reads in Hi-C data. Generative model to estimate probabilities for each bin-pair originating from a given origin. Reproducibility of contact matrices (stratum-adjusted correlation), reproducibility and number of significant interactions is improved. Novel interactions. Enrichment of TAD boundaries in LINE and SINE repetitive elements. Multi-mapping not sensitive to trimming. Read filtering strategy (Figure 1, supplementary figures are very visual). https://github.com/keleslab/mHiC
    - Zheng, Ye, Ferhat Ay, and Sunduz Keles. “Generative Modeling of Multi-Mapping Reads with MHi-C Advances Analysis of High Throughput Genome-Wide Conformation Capture Studies,” October 3, 2018. https://doi.org/10.1101/301705.

- `my5C`- web-based tools, well-documented analysis and visualization of 5S data, http://my5c.umassmed.edu/

- `nf-core-hic` - Analysis of Chromosome Conformation Capture data (Hi-C and more), Nextflow pipeline. https://github.com/nservant/nf-core-hic

- `TADbit` - Python-based pipeline, from iterative mapping, filtering, normalization. Similarity metrics: distance-centric Spearman, first principal eigenvector. TAD detection. TAD boundaries alignment, within 20kb. 3D modeling. Supplementary material - key functions, TAD detection algorithm, boundary comparison. https://github.com/3DGenomes/tadbit
    - Serra, François, Davide Baù, Mike Goodstadt, David Castillo, Guillaume J. Filion, and Marc A. Marti-Renom. “Automatic Analysis and 3D-Modelling of Hi-C Data Using TADbit Reveals Structural Features of the Fly Chromatin Colors.” PLoS Computational Biology 13, no. 7 (July 2017): e1005665. https://doi.org/10.1371/journal.pcbi.1005665.

### Mirnylab tools

- `cooler` file format for storing Hi-C matrices, sparse, hierarchical, multi-resolution. `cooler` Python package for data loading, aggregation, merging, normalization (balancing), viewing, exporting data. Together with "pairs" text-based format (https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md), and hic, cooler is accepted by the 4D Nucleome consortium DAC.https://github.com/mirnylab/cooler,https://cooler.readthedocs.io/en/latest/
    - Abdennur, Nezar, and Leonid Mirny. “Cooler: Scalable Storage for Hi-C Data and Other Genomically-Labeled Arrays.” BioRxiv, February 22, 2019. https://doi.org/10.1101/557660.

- `distiller-nf` - Java modular Hi-C mapping pipeline for reproducible data analysis, nextflow pipeline. Alignment, filtering, aggregating Hi-C matrices. https://github.com/mirnylab/distiller-nf

- `hiclib` - Python tools to qc, map, normalize, filter and analyze Hi-C data, https://bitbucket.org/mirnylab/hiclib

- `hic2cool` - Lightweight converter between hic and cool contact matrices. https://github.com/4dn-dcic/hic2cool

- `ICE` - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. 
    - Imakaev, Maxim, Geoffrey Fudenberg, Rachel Patton McCord, Natalia Naumova, Anton Goloborodko, Bryan R. Lajoie, Job Dekker, and Leonid A. Mirny. “Iterative Correction of Hi-C Data Reveals Hallmarks of Chromosome Organization.” Nature Methods 9, no. 10 (October 2012): 999–1003. https://doi.org/10.1038/nmeth.2148. - ICE - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. Assumption - all loci should have equal visibility. Deconvolution into eigenvectors/values. hiclib https://bitbucket.org/mirnylab/hiclib. Good description of the algorithm by Lior Pachter https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/

- `pairtools` - tools for low-level processing of mapped Hi-C paired reads. https://github.com/mirnylab/pairtools. Documentation, https://pairtools.readthedocs.io/en/latest/index.html

### Capture-C

- `CCseqBasic` - Capture-C analysis pipeline by Hughes lab, https://github.com/Hughes-Genome-Group/CCseqBasicS
- `CapSequm` - oligo design tool by Hughes lab, http://apps.molbiol.ox.ac.uk/CaptureC/cgi-bin/CapSequm.cgi

- `GOPHER` - probe design for Capture Hi-C. All, or selected, promoters, or around GWAS hits. Two other tools, CapSequm and HiCapTools. https://github.com/TheJacksonLaboratory/Gopher
    - Hansen, Peter, Salaheddine Ali, Hannah Blau, Daniel Danis, Jochen Hecht, Uwe Kornak, Darío G. Lupiáñez, Stefan Mundlos, Robin Steinhaus, and Peter N. Robinson. “GOPHER: Generator Of Probes for Capture Hi-C Experiments at High Resolution.” BMC Genomics 20, no. 1 (December 2019). https://doi.org/10.1186/s12864-018-5376-4.

- `capC-MAP` - Capture-C analysis pipeline. Python and C++, run through a configuration file. Outputs bedGraph. Compared with HiC-Pro, better detects PCR duplicates, identifies more interactions. Normalization tuned for Capture-C data.https://github.com/cbrackley/capC-MAP, https://capc-map.readthedocs.io/
    - Buckle, Adam, Nick Gilbert, Davide Marenduzzo, and Chris A Brackley. “CapC-MAP: A Software Package for Analysis of Capture-C Data.” Preprint. Genomics, October 30, 2018. https://doi.org/10.1101/456160.

### HiChIP

- `CID` - Chromatin Interaction Discovery, call chromatin interactions from ChIA-PET. Outperforms ChIA-PET2, MANGO pipelines, call more peaks than HICCUPS, hichipper. Java implementation, https://groups.csail.mit.edu/cgs/gem/cid/
    - Guo, Yuchun, Konstantin Krismer, Michael Closser, Hynek Wichterle, and David K Gifford. “High Resolution Discovery of Chromatin Interactions.” Nucleic Acids Research, February 14, 2019. https://doi.org/10.1093/nar/gkz051.

- `HiChIP-Peak` - HiChIP peak caller, focus on peaks at re-ligation sites. Peak filtering, then negative binomial model. Differential peak analysis similar to DiffBind. https://github.com/ChenfuShi/HiChIP_peaks
    - Shi, Chenfu, Magnus Rattray, and Gisela Orozco. “HiChIP-Peaks: A HiChIP Peak Calling Algorithm.” Preprint. Bioinformatics, June 27, 2019. https://doi.org/10.1101/682781.

### 4C

- `4Cseqpipe` processing pipeline and a genome-wide 4C primer database,  http://compgenomics.weizmann.ac.il/tanay/?page_id=367/
    - Werken, Harmen J. G. van de, Gilad Landan, Sjoerd J. B. Holwerda, Michael Hoichman, Petra Klous, Ran Chachik, Erik Splinter, et al. “Robust 4C-Seq Data Analysis to Screen for Regulatory DNA Interactions.” Nature Methods 9, no. 10 (October 2012): 969–72. https://doi.org/10.1038/nmeth.2173. - 4C technology paper. Two different 4bp cutters to increase resolution. Investigation of beta-globin locus, interchromosomal interactions.

### CUT&RUN

- CUT&RUN technology, chromatin profiling strategy, antibody-targeted controlled cleavage by micrococcal nuclease. Cost-efficient, low input requirements, easier.
    - Skene, Peter J, and Steven Henikoff. “An Efficient Targeted Nuclease Strategy for High-Resolution Mapping of DNA Binding Sites.” Genes and Chromosomes, n.d., 35. https://elifesciences.org/articles/21856

- `CUT&RUNTools` - a pipeline to fully process CUT&RUN data and identify protein binding and genomic footprinting from antibody-targeted primary cleavage data. Implemented in R, Python, Bach, runs under the SLURM job submission. At the core, creates a cut matrix of from enzyme cleavage data. Compared with Atactk and Centipede. https://bitbucket.org/qzhudfci/cutruntools/src/master/
    - Zhu, Qian. “CUT&RUNTools: A Flexible Pipeline for CUT&RUN Processing and Footprint Analysis,” 2019, 12. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4


## Resolution improvement


- `DeepHiC` - a generative adversarial network (GAN) for enhancing Hi-C data. Does not change the bin size, enhances the content of Hi-C data. Reconstructs the content from \~1% of the original data. Outperforms BoostHiC, HiCPlus, HiCNN. Online tool: http://sysomics.com/deephic/, code: https://github.com/omegahh/DeepHiC
    - Hong, Hao, Shuai Jiang, Hao Li, Cheng Quan, Chenghui Zhao, Ruijiang Li, Wanying Li, et al. “DeepHiC: A Generative Adversarial Network for Enhancing Hi-C Data Resolution.” Preprint. Bioinformatics, July 29, 2019. https://doi.org/10.1101/718148.

- `HiCNN` - a computational method for resolution enhancement. A modification of HiCPlus approach, using very deep (54 layers, five types of layers) convolutional neural network. A Hi-C matrix of regular resolution is transformed into high-resolution but very sparse matrix, HiCNN predicts the missing values. Pearson and MSE evaluation metrics, overlap of Fit-Hi-C-detected significant interactions - performs similar or slightly better than HiCPlus. PyTorch implementation. http://dna.cs.miami.edu/HiCNN/
    - Liu, Tong, and Zheng Wang. “HiCNN: A Very Deep Convolutional Neural Network to Better Enhance the Resolution of Hi-C Data.” Edited by John Hancock. Bioinformatics, April 9, 2019, btz251. https://doi.org/10.1093/bioinformatics/btz251.

- `Boost-HiC` - infer fine-resolution contact frequencies in Hi-C data, performs well even on 0.1% of the raw data. TAD boundaries remain. Better than HiCPlus. Can be used for differential analysis (comparison) of two Hi-C maps. https://github.com/LeopoldC/Boost-HiC
    - Carron, Leopold, Jean-baptiste Morlot, Vincent Matthys, Annick Lesne, and Julien Mozziconacci. “Boost-HiC : Computational Enhancement of Long-Range Contacts in Chromosomal Contact Maps,” November 18, 2018. https://doi.org/10.1101/471607.

- `HIFI` - Hi-C Interaction Frequency Inference for restriction fragment-resolution analysis of Hi-C data. Sparsity is resolved by using dependencies between neighboring restriction fragments, with Markov Random Fields performing the best. Better resolves TADs and sub-TADs, significant interactions. CTCF, RAD21, SMC3, ZNF143 are enriched around TAD boundaries. Matrices normalized for fragment-specific biases. https://github.com/BlanchetteLab/HIFI
    - Cameron, Christopher JF, Josée Dostie, and Mathieu Blanchette. “Estimating DNA-DNA Interaction Frequency from Hi-C Data at Restriction-Fragment Resolution.” Preprint. Bioinformatics, July 25, 2018. https://doi.org/10.1101/377523.

- `HiCPlus` - increasing resolution of Hi-C data using convolutional neural network. Basically, smoothing parts of Hi-C image, then binning into smaller parts. Performs better than bilinear/biqubic smoothing. https://github.com/zhangyan32/HiCPlus
    - Zhang, Yan, Lin An, Ming Hu, Jijun Tang, and Feng Yue. “HiCPlus: Resolution Enhancement of Hi-C Interaction Heatmap,” March 1, 2017. https://doi.org/10.1038/s41467-018-03113-2.


## Normalization

- `Binless` - a resolution-agnostic normalization method that adapts to the quality and quantity of available data, to detect significant interactions and differences. Negative binomial count regression framework, adapted for ICE normalization. Fused lasso to smooth neighboring signal. TADbit for data processing, details of read filtering. https://github.com/3DGenomes/binless
    - Spill, Yannick G., David Castillo, Enrique Vidal, and Marc A. Marti-Renom. “Binless Normalization of Hi-C Data Provides Significant Interaction and Difference Detection Independent of Resolution.” Nature Communications 10, no. 1 (26 2019): 1938. https://doi.org/10.1038/s41467-019-09907-2.

- `HiCNorm` - removing known biases in Hi-C data (GC content, mappability, fragment length) via Poisson regression, http://www.people.fas.harvard.edu/~junliu/HiCNorm/
    - Hu, Ming, Ke Deng, Siddarth Selvaraj, Zhaohui Qin, Bing Ren, and Jun S. Liu. “HiCNorm: Removing Biases in Hi-C Data via Poisson Regression.” Bioinformatics (Oxford, England) 28, no. 23 (December 1, 2012): 3131–33. https://doi.org/10.1093/bioinformatics/bts570. - Poisson normalization. Also tested negative binomial.

- `HiCorr` - visibility normalization using trans interactions only, to better emphasize promoter-enhancer interactions, combines advantages of implicit and explicit bias correction methods. https://github.com/JinLabBioinfo/HiCorr
    - Lu, Leina, Xiaoxiao Liu, Wei-Kai Huang, Paola Giusti-Rodriguez, Jian Cui, Shanshan Zhang, Wanying Xu, et al. “Robust Hi-C Chromatin Loop Maps in Human Neurogenesis and Brain Tissues at High-Resolution.” Preprint. Genomics, August 22, 2019. https://doi.org/10.1101/744540.

- `HiFive` - handling and normalization or pre-aligned Hi-C and 5C data, https://www.taylorlab.org/software/hifive/
    - Sauria, Michael EG, Jennifer E. Phillips-Cremins, Victor G. Corces, and James Taylor. “HiFive: A Tool Suite for Easy and Efficient HiC and 5C Data Analysis.” Genome Biology 16, no. 1 (December 2015). https://doi.org/10.1186/s13059-015-0806-y. - HiFive - post-processing of aligned Hi-C and 5C data, three normalization approaches: "Binning" - model-based Yaffe & Tanay's method, "Express" - matrix-balancing approach, "Probability" - multiplicative probability model. Judging normalization quality by correlation between matrices. 

- `HiTC` - The HiTC R package was developed to explore high-throughput 'C' data such as 5C or Hi-C. Dedicated R classes as well as standard methods for quality controls, normalization, visualization, and further analysis are also provided. https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization using loess regression on genomic distance, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. 


### CNV-aware normalization

- Hi-C data normalization considering CNVs. Extension of matrix-balancing algorithmto either retain the copy-number variation effect (LOIC) or remove them (CAIC). ICE itself can lead to misrepresentation of the contact probabilities between CNV regions. Estimating CNV directly from Hi-C data correcting for GC content, mappability, fragment length using Poisson regression. LOIC - the sum of contacts for a given genomic bin is proportional to CNV. CAIC - raw interaction counts are the product of a CNV bias matrix and the expected contact counts at a given genomic distance. Data, http://members.cbio.mines-paristech.fr/~nvaroquaux/normalization/, and `cancer-hic-norm` - Normalization of cancer Hi-C data, scripts for the manuscript. https://github.com/nservant/cancer-hic-norm. LOIC and CAIC methods are implemented in `iced` Python package, https://github.com/hiclib/iced
    - Servant, Nicolas, Nelle Varoquaux, Edith Heard, Emmanuel Barillot, and Jean-Philippe Vert. “Effective Normalization for Copy Number Variation in Hi-C Data.” BMC Bioinformatics 19, no. 1 (September 6, 2018): 313. https://doi.org/10.1186/s12859-018-2256-5.

- `HiCapp` - Iterative correction-based caICB method. Method to adjust for the copy number variants in Hi-C data. Loess-like idea - we converted the problem of removing the biases across chromosomes to the problem of minimizing the differences across count-distance curves of different chromosomes. Our method assumes equal representation of genomic locus pairs with similar genomic distances located on different chromosomes if there were no bias in the Hi-C maps. https://bitbucket.org/mthjwu/hicapp
    - Wu, Hua-Jun, and Franziska Michor. “A Computational Strategy to Adjust for Copy Number in Tumor Hi-C Data.” Bioinformatics (Oxford, England) 32, no. 24 (December 15, 2016): 3695–3701. https://doi.org/10.1093/bioinformatics/btw540.

- `OneD` - CNV bias-correction method, addresses the problem of partial aneuploidy. Bin-centric counts are modeled using negative binomial distribution, and its parameters are estimated using splines. A hidden Markov model is fit to infer copy number for each bin. Each Hi-C matrix entry is corrected by dividing its value by square root of the product of CNVs for the corresponding bins. Reproducibility score (eigenvector decomposition and comparison) to measure improvement in the similarity between replicated Hi-C data. https://github.com/qenvio/dryhic
    - Vidal, Enrique, François le Dily, Javier Quilez, Ralph Stadhouders, Yasmina Cuartero, Thomas Graf, Marc A Marti-Renom, Miguel Beato, and Guillaume J Filion. “OneD: Increasing Reproducibility of Hi-C Samples with Abnormal Karyotypes.” Nucleic Acids Research, January 31, 2018. https://doi.org/10.1093/nar/gky064.
    

## Reproducibility

- `“IDR2D` - Irreproducible Discovery Rate that identifies replicable interactions in ChIP-PET, HiChIP, and Hi-C data. Includes the original 1D IDR version (https://github.com/nboley/idr). Resolves multiple pairwise interactions.  https://github.com/gifford-lab/idr2d
    - Krismer, Konstantin, Yuchun Guo, and David K Gifford. “IDR2D Identifies Reproducible Genomic Interactions.” Preprint. Bioinformatics, July 3, 2019. https://doi.org/10.1101/691295.

- `3DChromatin_ReplicateQC` - Comparison of four Hi-C reproducibility assessment tools, `HiCRep`, `GenomeDISCO`, `HiC-Spector`, `QuASAR-Rep`. Tested the effects of noise, sparsity, resolution. Spearman doesn't work well. All tools performed similarly, worsening expectedly. QuASAR has QC tool measuring the level of noise. https://github.com/kundajelab/3DChromatin_ReplicateQC
    - Yardimci, Galip, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data,” September 14, 2017. doi:10.1101/188755. 

- `QuASAR` - Hi-C quality and reproducibility measure using spatial consistency between local and regional signals. Finds the maximum useful resolution by comparing quality and replicate scores of replicates. Part of HiFive pipeline,https://github.com/bxlab/hifive
    - Sauria, Michael EG, and James Taylor. “QuASAR: Quality Assessment of Spatial Arrangement Reproducibility in Hi-C Data.” BioRxiv, November 14, 2017. https://doi.org/10.1101/204438.

- `HiC-Spector` - reproducibility metric to quantify the similarity between contact maps using spectral decomposition. Decomposing Laplacian matrices and sum the Euclidean distance between eigenvectors. https://github.com/gersteinlab/HiC-spector
    - Yan, Koon-Kiu, Galip Gürkan Yardimci, Chengfei Yan, William S. Noble, and Mark Gerstein. “HiC-Spector: A Matrix Library for Spectral and Reproducibility Analysis of Hi-C Contact Maps.” Bioinformatics (Oxford, England) 33, no. 14 (July 15, 2017): 2199–2201. https://doi.org/10.1093/bioinformatics/btx152.

- `localtadsim` - Analysis of TAD similarity using variation of information (VI) metric as a local distance measure. 23 human Hi-C datasets, Hi-C Pro processed into 100kb matrices, Armatus to call TADs. Defining structurally similar and variable regions. Comparison with previous studies of genomic similarity. Cancer-normal comparison - regions containing pan-cancer genes are structurally conserved in normal-normal pairs, not in cancer-cancer. https://github.com/Kingsford-Group/localtadsim
    - Sauerwald, Natalie, and Carl Kingsford. “Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types.” Bioinformatics (Oxford, England) 34, no. 13 (July 1, 2018): i475–83. https://doi.org/10.1093/bioinformatics/bty265.


## Significant interaction (peak) callers

- `CHiCAGO` is a Capture Hi-C data processing method that filters out contacts that are expected by chance given the linear proximity of the interacting fragments on the genome and takes into account the asymmetric biases introduced by the capture step used in the Capture Hi-C approach. Two-component background model (Delaporte distribution) - Brownian motion (Neg. Binom.) and technical noise (Poisson). Account for distance. https://bioconductor.org/packages/release/bioc/html/Chicago.html
    - Cairns, Jonathan, Paula Freire-Pritchett, Steven W. Wingett, Csilla Várnai, Andrew Dimond, Vincent Plagnol, Daniel Zerbino, et al. “CHiCAGO: Robust Detection of DNA Looping Interactions in Capture Hi-C Data.” Genome Biology 17, no. 1 (2016): 127. https://doi.org/10.1186/s13059-016-0992-2.

- `ChiCMaxima` - a pipeline for detection and visualization of chromatin loops in Capture Hi-C data. Loess smoothing combined with a background model to detect significant interactions Comparison with GOTHiC and CHiCAGO. https://github.com/yousra291987/ChiCMaxima
    - Ben Zouari, Yousra, Anne M Molitor, Natalia Sikorska, Vera Pancaldi, and Tom Sexton. “ChiCMaxima: A Robust and Simple Pipeline for Detection and Visualization of Chromatin Looping in Capture Hi-C,” October 16, 2018. https://doi.org/10.1101/445023.

- `cLoops` - DBSCAN-based algorithm for the detection of chromatin loops in ChIA-PET, Hi-C, HiChIP, Trac-looping data. Local permutation-based estimation of statistical significance, several tests for enrichment over background. Outperforms diffHiC, Fit-Hi-C, GOTHiC, HiCCUPS, HOMER. https://github.com/YaqiangCao/cLoops
    - Cao, Yaqiang, Xingwei Chen, Daosheng Ai, Zhaoxiong Chen, Guoyu Chen, Joseph McDermott, Yi Huang, and Jing-Dong J. Han. “Accurate Loop Calling for 3D Genomic Data with CLoops,” November 8, 2018. https://doi.org/10.1101/465849.

- `coolpup.py` - Pile-up (aggregation, averaging) analysis of Hi-C data (.cool format) for visualizing and identifying chromatin loops from several sparse datasets, e.g., single-cell. Visualization using plotpup.py script. Scripts for paper: https://github.com/Phlya/coolpuppy_paper/tree/master/Nagano, tool: https://github.com/Phlya/coolpuppy
    - Flyamer, Ilya M., Robert S. Illingworth, and Wendy A. Bickmore. “Coolpup.Py - a Versatile Tool to Perform Pile-up Analysis of Hi-C Data.” BioRxiv, January 1, 2019, 586537. https://doi.org/10.1101/586537.

- `FastHiC` - hidden Markov random field (HMRF)-based peak caller, fast and well performing. https://yunliweb.its.unc.edu/fasthic/
    - Xu, Zheng, Guosheng Zhang, Cong Wu, Yun Li, and Ming Hu. “FastHiC: A Fast and Accurate Algorithm to Detect Long-Range Chromosomal Interactions from Hi-C Data.” Bioinformatics (Oxford, England) 32, no. 17 (01 2016): 2692–95. https://doi.org/10.1093/bioinformatics/btw240.

- `FIREcaller` - an R package to detect frequently interacting regions (FIREs, <200Kb interactions). Within-sample (HiCNormCis) and cross-sample (quantile) normalization, converting FIRE counts to Z-scores, taking significant ones. Schmitt data https://yunliweb.its.unc.edu/FIREcaller/
    - Crowley, Cheynna, Yuchen Yang, Yunjiang Qiu, Benxia Hu, Hyejung Won, Bing Ren, Ming Hu, and Yun Li. “FIREcaller: An R Package for Detecting Frequently Interacting Regions from Hi-C Data.” Preprint. Bioinformatics, April 29, 2019. https://doi.org/10.1101/619288.

- `Fit-Hi-C` - Python tool for detection of significant chromatin interactions, https://noble.gs.washington.edu/proj/fit-hi-c/
    - Ay, Ferhat, Timothy L. Bailey, and William Stafford Noble. “Statistical Confidence Estimation for Hi-C Data Reveals Regulatory Chromatin Contacts.” Genome Research 24, no. 6 (June 2014): 999–1011. https://doi.org/10.1101/gr.160374.113. - Fit-Hi-C method, Splines to model distance dependence. Model mid-range interaction frequencies, decay with distance. Biases, methods for normalization. Two-step splines - use all dots for first fit, identify and remove outliers, second fit without outliers. Markers of boundaries - insulators, heterochromatin, pluripotent factors. CNVs are enriched in chromatin boundaries. Replication timing data how-to http://www.replicationdomain.com/. Validation Hi-C data. http://chromosome.sdsc.edu/mouse/hi-c/download.html

- `FitHiChIP` - significant peak caller in HiChIP and PLAC-seq data. Accounts for assay-specific biases as well as for the distance effect. 3D differential loops detection. Methods. https://github.com/ay-lab/FitHiChIP
    - Bhattacharyya, Sourya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay. “FitHiChIP: Identification of Significant Chromatin Contacts from HiChIP Data,” September 10, 2018. https://doi.org/10.1101/412833.

- `GoTHIC` - R package for peak calling in individual HiC datasets, while accounting for noise. https://www.bioconductor.org/packages/release/bioc/html/GOTHiC.html
    - Mifsud, Borbala, Inigo Martincorena, Elodie Darbo, Robert Sugar, Stefan Schoenfelder, Peter Fraser, and Nicholas M. Luscombe. “GOTHiC, a Probabilistic Model to Resolve Complex Biases and to Identify Real Interactions in Hi-C Data.” Edited by Mark Isalan. PLOS ONE 12, no. 4 (April 5, 2017): e0174744. https://doi.org/10.1371/journal.pone.0174744. - The GOTHiC (genome organisation through HiC) algorithm uses a simple binomial distribution model to simultaneously remove coveralge-associated biases in Hi-C data and detect significant interactions by assuming that the global background interaction frequency of two loci. Use of the Benjamini–Hochberg multiple-testing correction to control for the false discovery rate. 

- `HiCPeaks` - Python CPU-based implementation for BH-FDR and HICCUPS, two peak calling algorithms for Hi-C data, proposed by Rao et al 2014. Text-to-cooler Hi-C data converter, two scripts to call peaks, and one for visualization (creation of a .png file)

- `HOMER` - Perl scripts for normalization, visualization, significant interaction detection, motif discovery. Does not correct for bias. http://homer.ucsd.edu/homer/interactions/

- `HiCapTools` - Software package that can design sequence capture probes for targeted chromosome capture applications and analyse sequencing output to detect proximities involving targeted fragments. Two probes are designed for each feature while avoiding repeat elements and non-unique regions. The data analysis suite processes alignment files to report genomic proximities for each feature at restriction fragment level and is isoform-aware for gene features. Statistical significance of contact frequencies is evaluated using an empirically derived background distribution. https://github.com/sahlenlab/HiCapTools
    - Anandashankar Anil, Rapolas Spalinskas, Örjan Åkerborg, Pelin Sahlén; HiCapTools: a software suite for probe design and proximity detection for targeted chromosome conformation capture applications, Bioinformatics, Volume 34, Issue 4, 15 February 2018, Pages 675–677, https://doi.org/10.1093/bioinformatics/btx625

- `HMRFBayesHiC` - a hidden Markov random field-based Bayesian peak caller to identify long-range chromatin interactions from Hi-C data. Borrowing information from neighboring loci. Previous peak calling methods, Fit-Hi-C. Interactions between enhancers and promoters as a benchmark. https://yunliweb.its.unc.edu/HMRFBayesHiC/
    - Xu, Zheng, Guosheng Zhang, Fulai Jin, Mengjie Chen, Terrence S. Furey, Patrick F. Sullivan, Zhaohui Qin, Ming Hu, and Yun Li. “A Hidden Markov Random Field-Based Bayesian Method for the Detection of Long-Range Chromosomal Interactions in Hi-C Data.” Bioinformatics (Oxford, England) 32, no. 5 (01 2016): 650–56. https://doi.org/10.1093/bioinformatics/btv650.

- `X-SCNN` - prediction of significant Hi-C interactions at highly improved resolution using TFBSs, histone marks, DNAse data (WIG format). A Siamese Convolutional Neural Network (SCNN) - two subnetworks with shared parameters predicting true interactions. HiCCUPS calls as true interactions, the same number of no interactions (balanced dataset). Keras with TensorFlow backend. https://github.com/ernstlab/X-SCNN.
    - Jaroszewicz, Artur, and Jason Ernst. “An Integrative Approach for Fine-Mapping Chromatin Interactions.” Preprint. Bioinformatics, April 11, 2019. https://doi.org/10.1101/605576.


## Differential interactions

- Check https://bitbucket.org/mforcato/hictoolscompare.git, they have tools for TAD comparison, and simulated matrices.

- TADbit paper [@Serra:2017aa]. Supplementary material http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005665#sec020 + key functions, TAD detection algorithm, border comparison. https://github.com/3DGenomes/tadbit

- `AP` - aggregation preference - parameter, to quantify TAD heterogeneity. Call significant interactions within a TAD, cluster with DBSCAN, calculate weighted interaction density within each cluster, average. AP measures are reproducible. Comparison of TADs in Gm12878 and IMR90 - stable TADs change their aggregation preference, these changes correlate with LINEs, Lamin B1 signal. Can detect structural changes (block split) in TADs. https://github.com/XiaoTaoWang/TADLib
    - Wang, X.-T., Dong, P.-F., Zhang, H.-Y., and Peng, C. (2015). Structural heterogeneity and functional diversity of topologically associating domains in mammalian genomes. Nucleic Acids Research 43, 7237–7246.
    
- `Chicdiff` - differential interaction detection in Capture Hi-C data. Signal normalization based on CHiCAGO framework, differential testing using DESeq2. Accounting for distance effect by the Independent Hypothesis Testing (IHW) method to learn p-value weights based on distance to maximize the number of rejected null hypotheses. https://github.com/RegulatoryGenomicsGroup/chicdiff
    - Cairns, Jonathan, William R. Orchard, Valeriya Malysheva, and Mikhail Spivakov. “Chicdiff: A Computational Pipeline for Detecting Differential Chromosomal Interactions in Capture Hi-C Data.” BioRxiv, January 1, 2019, 526269. https://doi.org/10.1101/526269.

- `diffHiC` - Differential contacts using the full pipeline for Hi-C data. Explanation of the technology, binning. MA normalization, edgeR-based. Comparison with HOMER. https://bioconductor.org/packages/release/bioc/html/diffHic.html
    - Lun, Aaron T. L., and Gordon K. Smyth. “DiffHic: A Bioconductor Package to Detect Differential Genomic Interactions in Hi-C Data.” BMC Bioinformatics 16 (2015): 258. https://doi.org/10.1186/s12859-015-0683-0.

- `diffloop` - Differential analysis of chromatin loops (ChIA-PET). edgeR framework. https://bioconductor.org/packages/release/bioc/html/diffloop.html
    - Lareau, Caleb A., and Martin J. Aryee. “Diffloop: A Computational Framework for Identifying and Analyzing Differential DNA Loops from Sequencing Data.” Bioinformatics (Oxford, England), September 29, 2017. https://doi.org/10.1093/bioinformatics/btx623.

- DiffTAD - differential contact frequency in TADs between two conditions. Two - permutation-based comparing observed vs. expected median interactions, and parametric test considering the sign of the differences within TADs. Both tests account for distance stratum. https://bitbucket.org/rzaborowski/differential-analysis
    - Zaborowski, Rafal, and Bartek Wilczynski. “DiffTAD: Detecting Differential Contact Frequency in Topologically Associating Domains Hi-C Experiments between Conditions.” BioRxiv, January 1, 2016, 093625. https://doi.org/10.1101/093625.

- `FIND` - differential chromatin interaction detection comparing the local spatial dependency between interacting loci. Previous strategies - simple fold-change comparisons, binomial model (HOMER), count-based (edgeR). FIND exploits a spatial Poisson process model to detect differential chromatin interactions that show both a significant change in their interaction frequency and the interaction frequency of their adjacent bins. "Variogram" concept. For each point, compare densities between conditions using Fisher's test. Explored various multiple correction testing methods, used r^th ordered p-values (rOP) method. Benchmarking against edgeR in simulated settings - FIND outperforms at shorter distances, edgeR has more false positives at longer distances. Real Hi-C data normalized using KR and MA normalizations. R paclage https://bitbucket.org/nadhir/find/downloads/
    - Mohamed Nadhir, Djekidel, Yang Chen, and Michael Q. Zhang. “FIND: DifFerential Chromatin INteractions Detection Using a Spatial Poisson Process.” Genome Research, February 12, 2018. https://doi.org/10.1101/gr.212241.116.

- `HiCcompare` - joint normalization of two Hi-C datasets using loess regression through an MD plot (minus-distance). Data-driven normalization accounting for the between-dataset biases. Per-distance permutation testing of significant interactions. http://bioconductor.org/packages/release/bioc/html/HiCcompare.html
    - Stansfield, John C., Kellen G. Cresswell, Vladimir I. Vladimirov, and Mikhail G. Dozmorov. “HiCcompare: An R-Package for Joint Normalization and Comparison of HI-C Datasets.” BMC Bioinformatics 19, no. 1 (December 2018). https://doi.org/10.1186/s12859-018-2288-x.

- `multiHiCcompare` - joint normalization of multiple Hi-C datasets using cyclic loess regression through pairs of MD plots (minus-distance). Data-driven normalization accounting for the between-dataset biases. Per-distance edgeR-based testing of significant interactions. http://bioconductor.org/packages/release/bioc/html/multiHiCcompare.html
    - Stansfield, John C, Kellen G Cresswell, and Mikhail G Dozmorov. “MultiHiCcompare: Joint Normalization and Comparative Analysis of Complex Hi-C Experiments.” Edited by Inanc Birol. Bioinformatics, January 22, 2019. https://doi.org/10.1093/bioinformatics/btz048.

- `Selfish` - comparative analysis of replicate Hi-C experiments via a self-similarity measure - local similarity borrowed from image comparison. Check reproducibility, detect differential interactions. Boolean representation of contact matrices for reproducibility quantification. Deconvoluting local interactions with a Gaussian filter (putting a Gaussian bell around a pixel), then comparing derivatives between contact maps for each radius. Simulated (Zhou method) and real comparison with FIND - better performance, especially on low fold-changes. Stronger enrichment of relevant epigenomic features. Matlab implementation https://github.com/ucrbioinfo/Selfish
    - Roayaei Ardakany, Abbas, Ferhat Ay, and Stefano Lonardi. “Selfish: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure.” BioRxiv, January 1, 2019, 540708. https://doi.org/10.1101/540708.



## TAD callers

- [Brief description of 22 TAD calling methods](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1596-9/MediaObjects/13059_2018_1596_MOESM1_ESM.pdf). Source: [Zufferey et al., “Comparison of Computational Methods for the Identification of Topologically Associating Domains.”](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1596-9#Bib1)

- `TADpole` - hierarchical TAD boundary caller. Preprocessing by filtering sparse rows, transforming the matrix into its Pearson correlation coefficient matrix, running PCA on it and retaining 200 PCs, transforming into a Euclidean distance matrix, clustering using the Constrained Incremental Sums of Squares clustering (rioja::chclust(, coniss)), estimating significance, Calinski-Harabasz index to estimate the optimal number of clusters (chromatin subdivisions). Benchmarking using Zufferey 2018 datasets, mouse limb bud development with genomic inversions from Kraft 2019. Resolution, normalization, sequencing depth. Metrics: the Overlap Score, the Measure of Concordance, all from Zufferey 2018. Enrichment in epigenomic marks. DiffT metric for differential analysis (on binarized TAD/non-TAD matrices). Compared with 22 TAD callers, including hierarchical (CaTCH, GMAP, Matryoshka, PSYCHIC). https://github.com/3DGenomes/TADpole
    - Soler-Vila, Paula, Pol Cuscó Pons, Irene Farabella, Marco Di Stefano, and Marc A. Marti-Renom. “Hierarchical Chromatin Organization Detected by TADpole.” Preprint. Bioinformatics, July 11, 2019. https://doi.org/10.1101/698720.

- `HiCDB` - TAD boundary detection using local relative insulation (LRI) metric, improved stability, less parameter tuning, cross-resolution, differential boundary detection, lower computations, visualization. Review of previous methods, directionality index, insulation score. Math of LRI. GSEA-like enrichment in genome annotations (CTCF). Differential boundary detection using intersection of extended boundaries. Compared with Armatus, DI, HiCseg, IC-finder, Insulation, TopDom on 40kb datasets. Accurately detects smaller-scale boundaries. Differential TADs are enriched in cell type-specific genes. https://github.com/ChenFengling/RHiCDB
    - Chen, Fengling, Guipeng Li, Michael Q. Zhang, and Yang Chen. “HiCDB: A Sensitive and Robust Method for Detecting Contact Domain Boundaries.” Nucleic Acids Research 46, no. 21 (November 30, 2018): 11239–50. https://doi.org/10.1093/nar/gky789.

- `OnTAD` - hierarchical TAD caller, Optimal Nested TAD caller. Sliding window, adaptive local minimum search algorithm, similar to TOPDOM. https://github.com/anlin00007/OnTAD
    - An, Lin, Tao Yang, Jiahao Yang, Johannes Nuebler, Qunhua Li, and Yu Zhang. “Hierarchical Domain Structure Reveals the Divergence of Activity among TADs and Boundaries,” July 3, 2018. https://doi.org/10.1101/361147. - Intro about TADs, Dixon's directionality index, Insulation score. Other hierarchical callers - TADtree, rGMAP, Arrowhead, 3D-Net, IC-Finder. Limitations of current callers - ad hoc thresholds, sensitivity to sequencing depth and mapping resolution, long running time and large memory usage, insufficient performance evaluation. Boundaries are asymmetric - some have more contacts with other boundaries, support for asymmetric loop extrusion model. Performance comparison with DomainCaller, rGMAP, Arrowhead, TADtree. Stronger enrichment of CTCF and two cohesin proteins RAD21 and SMC3. TAD-adjR^2 metric quantifying proportion of variance in the contact frequencies explained by TAD moundaries. Reproducibility of TAD boundaries - Jaccard index, tested at different sequencing depths and resolutions. Boundaries of hierarchical TADs are more active - more CTCF, epigenomic features, TFBSs expressed genes. Super-boundaries - shared by 5 or more TADs, highly active. Rao-Huntley 2014 Gm12878 data. Distance correction - subtracting the mean counts at each distance.

- `3D-NetMod` - hierarchical, nested, partially overlapping TAD detection using graph theory. Community detection method based on the maximization of network modularity, Louvain-like locally greedy algorithm, repeated several (20) times to avoid local maxima, then getting consensus. Tuning parameters are estimated over sequence search. Benchmarked against TADtree, directionality index, Arrowhead. ICE-normalized data brain data from Geschwind (human data) and Jiang (mouse data) studies. Computationally intensive. Python implementation https://bitbucket.org/creminslab/3dnetmod_method_v1.0_10_06_17
    - Norton, Heidi K., Daniel J. Emerson, Harvey Huang, Jesi Kim, Katelyn R. Titus, Shi Gu, Danielle S. Bassett, and Jennifer E. Phillips-Cremins. “Detecting Hierarchical Genome Folding with Network Modularity.” Nature Methods 15, no. 2 (February 2018): 119–22. https://doi.org/10.1038/nmeth.4560.

- `deDoc` - TAD detection minimizing structural entropy of the Hi-C graph (structural information theory). Detects optimal resolution (= minimal entropy). Pooled 10 single-cell Hi-C analysis. Intro about TADs, brief description of TAD callers, including hierarchical. Works best on raw, non-normalized data, highly robust to sparsity (0.1% of the original data sufficient). Compared with five TAD callers (Armatus, TADtree, Arrowhead, MrTADFinder, Domaincall (DI)), and a classical graph modularity detection algorithm. Enrichment in CTCF, housekeeping genes, H3K4me3, H4K20me1, H3K36me3. Other benchmarks - weighted similarity, number, length of TADs. Detects hierarchy over different passes. Java implementation (wont run on Mac)https://github.com/yinxc/structural-information-minimisation
    - Li, Angsheng, Xianchen Yin, Bingxiang Xu, Danyang Wang, Jimin Han, Yi Wei, Yun Deng, Ying Xiong, and Zhihua Zhang. “Decoding Topologically Associating Domains with Ultra-Low Resolution Hi-C Data by Graph Structural Entropy.” Nature Communications 9, no. 1 (15 2018): 3265. https://doi.org/10.1038/s41467-018-05691-7.

- `CaTCH` - identification of hierarchical TAD structure, https://github.com/zhanyinx/CaTCH_R
    - Zhan, Yinxiu, Luca Mariani, Iros Barozzi, Edda G. Schulz, Nils Blüthgen, Michael Stadler, Guido Tiana, and Luca Giorgetti. “Reciprocal Insulation Analysis of Hi-C Data Shows That TADs Represent a Functionally but Not Structurally Privileged Scale in the Hierarchical Folding of Chromosomes.” Genome Research 27, no. 3 (2017): 479–90. https://doi.org/10.1101/gr.212803.116. - CaTCH - identification of hierarchical TAD structure. Reciprocal insulation (RI) index. Benchmarked against Dixon's TADs (diTADs). CTCF enrichment as a benchmark, enrichment of TADs in differentially expressed genes. https://github.com/zhanyinx/CaTCH_R

- `HiTAD` - hierarchical TAD identification, different resolutions, correlation with chromosomal compartments, replication timing, gene expression. Adaptive directionality index approach. Data sources, methods for comparing TAD boundaries, reproducibility. H3K4me3 enriched and H3K4me1 depleted at boundaries. TAD boundaries (but not sub-TADs) separate replication timing, A/B compartments, gene expression. https://github.com/XiaoTaoWang/TADLib,  https://pypi.python.org/pypi/TADLib
    - Wang, Xiao-Tao, Wang Cui, and Cheng Peng. “HiTAD: Detecting the Structural and Functional Hierarchies of Topologically Associating Domains from Chromatin Interactions.” Nucleic Acids Research 45, no. 19 (November 2, 2017): e163. https://doi.org/10.1093/nar/gkx735.

- `IC-Finder` - Segmentations of HiC maps into hierarchical interaction compartments, http://membres-timc.imag.fr/Daniel.Jost/DJ-TIMC/Software.html
    - Noelle Haddad, Cedric Vaillant, Daniel Jost. "IC-Finder: inferring robustly the hierarchical organization of chromatin folding." Nucleic Acids Res. 2017 Jun 2; 45(10). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5449546/. TAD-identification procedure that is essentially identical to my idea, including inferring hierarchical structure of TADs based on hierarchical clustering. Article does not explicity discuss A/B compartments so it may still be novel for compartment detection. - TAD finding by combined hierarchical clustering. Correlation-based distance, weighted-mean linkage similarity. Variance criterion to define cluster boundaries. http://membres-timc.imag.fr/Daniel.Jost/DJ-TIMC/Software.html

- `ClusterTAD` - A clustering method for identifying topologically associated domains (TADs) from Hi-C data, https://github.com/BDM-Lab/ClusterTAD
    - Oluwadare, Oluwatosin, and Jianlin Cheng. “ClusterTAD: An Unsupervised Machine Learning Approach to Detecting Topologically Associated Domains of Chromosomes from Hi-C Data.” BMC Bioinformatics 18, no. 1 (November 14, 2017): 480. https://doi.org/10.1186/s12859-017-1931-2. https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1931-2 - ClusterTAD paper. Clustering to define TADs. Datasets: simulated Hi-C data with pre-defined TADs https://link.springer.com/article/10.1007%2Fs40484-015-0047-9, RenLab Hi-C and CTCF data http://chromosome.sdsc.edu/mouse/download.html. 

- `EAST` - Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps, https://github.com/ucrbioinfo/EAST
    - Abbas Roayaei Ardakany, Stefano Lonardi, and Marc Herbstritt, “Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps” (Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik GmbH, Wadern/Saarbruecken, Germany, 2017), https://doi.org/10.4230/LIPIcs.WABI.2017.22. - EAST: Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps. Haar-like features (rectangles on images) and a function that quantifies TAD properties: frequency within is high, outside - low, boundaries must be strong. Objective - finding a set of contigious non-overlapping domains maximizing the function. Restricted by maximum length of TADs Text boundaries for enrichment in CTCF, RNP PolII, H3K4me3, H3K27ac. https://github.com/ucrbioinfo/EAST

- `TADtree` - TADtree is an algorithm the identification of hierarchical topological domains in Hi-C data, http://compbio.cs.brown.edu/software/
    - Weinreb, Caleb, and Benjamin J. Raphael. “Identification of Hierarchical Chromatin Domains.” Bioinformatics (Oxford, England) 32, no. 11 (June 1, 2016): 1601–9. https://doi.org/10.1093/bioinformatics/btv485. - TADtree paper. Hierarchical (nested) TAD identification. Two ways of TAD definition: 1D and 2D. Normalization by distance. Enrichment over background. Deep statistics of the method. How to compare TADs (VI measure (`vi.dist` in https://cran.r-project.org/web/packages/mcclust/mcclust.pdf), Precision/recall using Dixon as the true set, Fig. 5: number of TADs, TAD size boxplots, Enrichment within 50kb of a TAD boundary - CTCF, PolII, H3K4me3, housekeeping genes - stronger enrichment the better). http://compbio.cs.brown.edu/software/

- `TopDom` - An efficient and Deterministic Method for identifying Topological Domains in Genomes, http://zhoulab.usc.edu/TopDom/
    - Shin, Hanjun, Yi Shi, Chao Dai, Harianto Tjong, Ke Gong, Frank Alber, and Xianghong Jasmine Zhou. “TopDom: An Efficient and Deterministic Method for Identifying Topological Domains in Genomes.” Nucleic Acids Research 44, no. 7 (April 20, 2016): e70. https://doi.org/10.1093/nar/gkv1505. - TopDom paper. Review of other methods. Method is based on general observation that within-TAD interactions are stronger than between-TAD. binSignal value as the average of nearby contact frequency, fitting a curve, finding local minima, test them for significance. Fast, takes linear time. Detects similar domains to HiCseq and Dixon's directionaliry index. Found expected enrichment in CTCF, histone marks. Housekeeping genes and overall gene density are close to TAD boundaries, differentially expressed genes are not. Figure 7 - how to detect common/unique boundaries using Jaccard-like statistics. http://zhoulab.usc.edu/TopDom/

- `Armatus` - TAD detection at different resolutions, https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus
    - Filippova, Darya, Rob Patro, Geet Duggal, and Carl Kingsford. “Identification of Alternative Topological Domains in Chromatin.” Algorithms for Molecular Biology 9, no. 1 (2014): 14. doi:10.1186/1748-7188-9-14.https://almob.biomedcentral.com/track/pdf/10.1186/1748-7188-9-14 - Dynamic programming method named “Armatus”. Methods - statistics, intuitive. Consider different resolutions, stability of the detected TADs, persistency across resolution. Only one parameter controlling the size of domains. Their TADs are different from Dixon's TADs. https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus

- `HiCseg` - TAD detection by maximization of likelihood based block-wise segmentation model, https://cran.r-project.org/web/packages/HiCseg/index.html
    - Lévy-Leduc, Celine, M. Delattre, T. Mary-Huard, and S. Robin. “Two-Dimensional Segmentation for Analyzing Hi-C Data.” Bioinformatics (Oxford, England) 30, no. 17 (September 1, 2014): i386-392. https://doi.org/10.1093/bioinformatics/btu443. - HiCseg paper. TAD detection by maximization of likelihood based block-wise segmentation model, HiCseg R package. 2D segmentation rephrased as 1D segmentation - not contours, but borders. Statistical framework, solved with dynamic programming. Dixon data as gold standard. Hausdorff distance to compare segmentation quality, https://en.wikipedia.org/wiki/Hausdorff_distance. Parameters (from TopDom paper): nb_change_max = 500, distrib = 'G' and model = 'Dplus'.

- `hickit` - TAD calling, phase imputation, 3D modeling and more for diploid single-cell Hi-C (Dip-C) and bulk Hi-C, https://github.com/lh3/hickit

- Insulation score: Giorgetti, Luca, Bryan R. Lajoie, Ava C. Carter, Mikael Attia, Ye Zhan, Jin Xu, Chong Jian Chen, et al. “Structural Organization of the Inactive X Chromosome in the Mouse.” Nature 535, no. 7613 (28 2016): 575–79. https://doi.org/10.1038/nature18589. - https://github.com/dekkerlab/cworld-dekker/tree/master/scripts/perl matrix2insulation.pl, Parameters: -is 480000 -ids 320000 -im iqrMean -nt 0 -ss 160000 -yb 1.5 -nt 0 -bmoe 0. 
    - Used in Yardımcı, Galip Gürkan, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data.” BioRxiv, January 1, 2018. https://doi.org/10.1101/188755. - TADs detected by insulation score are robust to resolution and noise


## Prediction of 3D features

- `3DEpiLoop` - prediction of 3D interactions from 1D epigenomic profiles using Random Forest trained on CTCF peaks (histone modifications are the most important predictors, and TFBSs). https://bitbucket.org/4dnucleome/3depiloop
    - Al Bkhetan, Ziad, and Dariusz Plewczynski. “Three-Dimensional Epigenome Statistical Model: Genome-Wide Chromatin Looping Prediction.” Scientific Reports 8, no. 1 (December 2018). https://doi.org/10.1038/s41598-018-23276-8.

- `SNIPER` - 3D subcompartment (A1, A2, B1, B2, B3) identification from low-coverage Hi-C datasets. A neural network based on a denoising autoencoder (9 layers) and multi-layer perceptron. Sigmoidal activation of inputs, ReLU, softmax on outputs. Dropout, binary cross-entropy. exp(-1/C) transformation of Hi-C matrices. Applied to Gm12878 and 8 additional cell types to compare subcompartment changes. Compared with Rao2014 annotations, outperforms Gaussian HMM and MEGABASE. https://github.com/ma-compbio/SNIPER
    - Xiong, Kyle, and Jian Ma. “Revealing Hi-C Subcompartments by Imputing High-Resolution Inter-Chromosomal Chromatin Interactions.” BioRxiv, January 1, 2018, 505503. https://doi.org/10.1101/505503.

- `TADBoundaryDectector` - TAD boundary prediction from sequence only using deep learning models. 12 architectures tested, with three convolutional and an LSTM layer performed best. Methods, Implementation in Keras-TensorFlow. Model evaluation using different criteria, 96% accuracy reported. Deep learning outperform feature-based models, among which Boosted trees, Random Forest, elastic net logistic regression are best performers. Data augmentation (aka feature engineering) by randomly shifting TAD boundary regions by som base pairs of length (0-100). Tested on Drozophila data. https://github.com/lincshunter/TADBoundaryDectector
    - Henderson, John, Vi Ly, Shawn Olichwier, Pranik Chainani, Yu Liu, and Benjamin Soibam. “Accurate Prediction of Boundaries of High Resolution Topologically Associated Domains (TADs) in Fruit Flies Using Deep Learning.” Nucleic Acids Research, May 3, 2019. https://doi.org/10.1093/nar/gkz315.


## SNP-oriented

- `iRegNet3D` - Integrated Regulatory Network 3D (iRegNet3D) is a high-resolution regulatory network comprised of interfaces of all known transcription factor (TF)-TF, TF-DNA interaction interfaces, as well as chromatin-chromatin interactions and topologically associating domain (TAD) information from different cell lines.  
    - Goal: SNP interpretation
    - Input: One or several SNPs, rsIDs or genomic coordinates.
    - Output: For one or two SNPs, on-screen information of their disease-related info, connection over TF-TF and chromatin interaction networks, and whether they interact in 3D and located within TADs. For multiple SNPs, same info downloadable as text files.

- `3DSNP` - A database linking noncoding SNPs to 3D interacting genes. http://cbportal.org/3dsnp/
    - Lu, Yiming, Cheng Quan, Hebing Chen, Xiaochen Bo, and Chenggang Zhang. “3DSNP: A Database for Linking Human Noncoding SNPs to Their Three-Dimensional Interacting Genes.” Nucleic Acids Research 45, no. D1 (January 4, 2017): D643–49. https://doi.org/10.1093/nar/gkw1022. - 3DSNP database integrating SNP epigenomic annotations with chromatin loops. Linear closest gene, 3D interacting gene, eQTL, 3D interacting SNP, chromatin states, TFBSs, conservation. For individual SNPs.

- HUGIn, tissue-specific Hi-C linear display of anchor position and around. Overlay gene expression and epigenomic data. Association of SNPs with genes based on Hi-C interactions. Tissue-specific. http://yunliweb.its.unc.edu/HUGIn/
    - Martin, Joshua S, Zheng Xu, Alex P Reiner, Karen L Mohlke, Patrick Sullivan, Bing Ren, Ming Hu, and Yun Li. “HUGIn: Hi-C Unifying Genomic Interrogator.” BioRxiv, 2017, 117531.


## CNV and Structural variant detection

- [hicpipe](https://github.com/ChenFengling/HiCpipe) and, alternatively, [HiCnorm](http://www.people.fas.harvard.edu/~junliu/HiCNorm/) normalization preserves CNVs in Hi-C data. From Zhang et al., “Local and Global Chromatin Interactions Are Altered by Large Genomic Deletions Associated with Human Brain Development.” https://doi.org/10.1038/s41467-018-07766-x

- `hic_breakfinder` - SV identification in Hi-C data. https://github.com/dixonlab/hic_breakfinder
    - Dixon, Jesse R., Jie Xu, Vishnu Dileep, Ye Zhan, Fan Song, Victoria T. Le, Galip Gürkan Yardımcı, et al. “Integrative Detection and Analysis of Structural Variation in Cancer Genomes.” Nature Genetics, September 10, 2018. https://doi.org/10.1038/s41588-018-0195-8. - Detection of structural variants (SV) by integrating optical mapping, Hi-C, and WGS. Custom pipeline using LUMPY, Delly, Control-FREEC software. New Hi-C data on 14 cancer cell lines and 21 previously published datasets. Integration of the detected SVs with genomic annotations, including replication timing. Supplementary data with SVs resolved by individual methods and integrative approaches.

- `HiCnv` - CNV, translocation calling from Hi-C data. CNV calling using HMM on per-restriction site quantified data and 1D-normalized accounting for low GC-content (<0.2), mappability (<0.5). Translocation calling on inter-chromosomal matrices, binned. CNV calling: https://github.com/ay-lab/HiCnv, Translocation calling: https://github.com/ay-lab/HiCtrans, Hi-C simulation: https://github.com/ay-lab/AveSim 
    - Chakraborty, Abhijit, and Ferhat Ay. “Identification of Copy Number Variations and Translocations in Cancer Cells from Hi-C Data.” Edited by Christina Curtis. Bioinformatics 34, no. 2 (January 15, 2018): 338–45. https://doi.org/10.1093/bioinformatics/btx664.

- `HiNT` - CNV and translocation detection from \~10-20% ambigious chimeric reads in Hi-C data. Three tools: HiNT-Pre - preprocessing of Hi-C data; HiNT-CNV and HiNT-TL - CNV and translocation detection, respectively (accept HiC-Pro output). Tested on K562 (cancer) and Gm12878 (normal) data. Removal of known biases using a GAM with Poisson function. Outperforms Delly, Meerkat, hic_breakfinder, HiCtrans. Relatively little overlap with CNVs from WGS (BIC-seq2). Gold-standard - FISH data from Dixon et al., “Integrative Detection and Analysis of Structural Variation in Cancer Genomes.” https://github.com/parklab/HiNT
    - Wang, Su, Soohyun Lee, Chong Chu, Dhawal Jain, Geoff Nelson, Jennifer M. Walsh, Burak H. Alver, and Peter J. Park. “HiNT: A Computational Method for Detecting Copy Number Variations and Translocations from Hi-C Data.” Preprint. Bioinformatics, June 3, 2019. https://doi.org/10.1101/657080.


## Visualization

- Hi-C data visualization review. Good introduction into the 3D genome organization, 115 key references. [Table 2. Hi-C visualization tools](https://dev.biologists.org/highwire/markup/1255595/expansion?width=1000&height=500&iframe=true&postprocessors=highwire_tables%2Chighwire_reclass%2Chighwire_figures%2Chighwire_math%2Chighwire_inline_linked_media%2Chighwire_embed)
    - Ing-Simmons, Elizabeth, and Juan M. Vaquerizas. “Visualising Three-Dimensional Genome Organisation in Two Dimensions.” Development 146, no. 19 (October 1, 2019): dev177162. https://doi.org/10.1242/dev.177162.

- `3D Genome Browser` - visualizing existing Hi-C and other chromatin conformation capture data. Alongside with genomic and epigenomic data. Own data can be submitted in BUTLR format. http://promoter.bx.psu.edu/hi-c/
    - Wang, Yanli, Bo Zhang, Lijun Zhang, Lin An, Jie Xu, Daofeng Li, Mayank NK Choudhary, et al. “The 3D Genome Browser: A Web-Based Browser for Visualizing 3D Genome Organization and Long-Range Chromatin Interactions.” BioRxiv, 2017, 112268.

- `CSynth` - 3D genome interactive modeling on GPU, and visualization. http://csynth.org/
    - Todd, Stephen, Peter Todd, Simon J McGowan, James R Hughes, Yasutaka Kakui, Frederic Fol Leymarie, William Latham, and Stephen Taylor. “CSynth: A Dynamic Modelling and Visualisation Tool for 3D Chromatin Structure.” BioRxiv, January 1, 2019, 499806. https://doi.org/10.1101/499806.

- `DNARchitect` - a Shiny App for visualizing genomic data (HiC, mRNA, ChIP, ATAC etc) in bed, bedgraph, and bedpe formats. Web version, http://shiny.immgen.org/DNARchitect/, GitHub, https://github.com/alosdiallo/DNA_Rchitect

- `GENOVA` - GENome Organisation Visual Analytics, an R package for rich visual analysis of Hi-C data. Input - HiC-Pro processed files, BED, text formats. Single or two experiment analysis. Integration of external annotations, A/B compartments, cis-/trans-interactions, TADs and loops, genes, insluation score heatmap, differences. https://github.com/robinweide/GENOVA

- `HiCExplorer` - set of programs to process, normalize, analyze and visualize Hi-C data, Python. https://hicexplorer.readthedocs.io/en/latest/, https://github.com/deeptools/HiCExplorer/. 
- ChoroGenome browser - genes vs. TAD boundaries. http://chorogenome.ie-freiburg.mpg.de/, and the underlying HiCBrowser https://github.com/deeptools/HiCBrowser/
    - Ramírez, Fidel, Vivek Bhardwaj, Laura Arrigoni, Kin Chung Lam, Björn A. Grüning, José Villaveces, Bianca Habermann, Asifa Akhtar, and Thomas Manke. “High-Resolution TADs Reveal DNA Sequences Underlying Genome Organization in Flies.” Nature Communications 9, no. 1 (December 2018). https://doi.org/10.1038/s41467-017-02525-w.
- `Galaxy HiCExplorer` - a web server for Hi-C data preprocessing, QC, visualization. Web interface,   https://hicexplorer.usegalaxy.eu/, Docker container, https://github.com/deeptools/docker-galaxy-hicexplorer
    - Wolff, Joachim, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach, Thomas Manke, Rolf Backofen, Fidel Ramírez, and Björn A. Grüning. “Galaxy HiCExplorer: A Web Server for Reproducible Hi-C Data Analysis, Quality Control and Visualization.” Nucleic Acids Research 46, no. W1 (July 2, 2018): W11–16. https://doi.org/10.1093/nar/gky504.


- `HiGlass` visualization server for Google maps-style navigation of Hi-C maps. Overlay genes, epigenomic tracks. http://higlass.io/, https://github.com/higlass/higlass, and many HiGlass-related developmend from the author, https://github.com/pkerpedjiev
    - Kerpedjiev, Peter, Nezar Abdennur, Fritz Lekschas, Chuck McCallum, Kasper Dinkla, Hendrik Strobelt, Jacob M Luber, et al. “HiGlass: Web-Based Visual Comparison And Exploration Of Genome Interaction Maps.” BioRxiv, 2017, 121889.
    - Python bindings to and Jupyter Notebook+Lab integration for the HiGlass viewer. https://github.com/higlass/higlass-python

- `HiPiler` - exploration and comparison of loops and domains as snippets-heatmaps of data. https://github.com/flekschas/hipiler
    - Lekschas, Fritz, Benjamin Bach, Peter Kerpedjiev, Nils Gehlenborg, and Hanspeter Pfister. “HiPiler: Visual Exploration of Large Genome Interaction Matrices with Interactive Small Multiples.” IEEE Transactions on Visualization and Computer Graphics 24, no. 1 (January 2018): 522–31. https://doi.org/10.1109/TVCG.2017.2745978. TechBlog: HiPiler simplifies chromatin structure analysis, http://blogs.nature.com/naturejobs/2017/09/11/techblog-hipiler-simplifies-chromatin-structure-analysis/

- `HiCPlotter` - Hi-C visualization tool, allows for integrating various data tracks. https://github.com/kcakdemir/HiCPlotter
    - Akdemir, Kadir Caner, and Lynda Chin. “HiCPlotter Integrates Genomic Data with Interaction Matrices.” Genome Biology 16 (2015): 198. https://doi.org/10.1186/s13059-015-0767-1.

- `HiTC` - plotting of rotated halves of chromosome interaction matrix, with annotations. https://bioconductor.org/packages/release/bioc/html/HiTC.html

- `NAT` - the 4D Nucleome Analysis Toolbox, for Hi-C data (text, cool format) normalization (ICE, Toeplitz, CNV-Toeplitz), TAD calling (Directionality index, Armatus, custom), karyotype abnormalities visualization on inter-chromosomal matrices, timecourse visualization. Matlab.https://github.com/laseaman/4D_Nucleome_Analysis_Toolbox
    - Seaman, Laura, and Indika Rajapakse. “4D Nucleome Analysis Toolbox: Analysis of Hi-C Data with Abnormal Karyotype and Time Series Capabilities.” Bioinformatics (Oxford, England) 34, no. 1 (01 2018): 104–6. https://doi.org/10.1093/bioinformatics/btx484.

- `pyGenomeTracks` - python module to plot beautiful and highly customizable genome browser tracks, https://github.com/deeptools/pyGenomeTracks

- TADKit - 3D Genome Browser. Main web site, http://sgt.cnag.cat/3dg/tadkit/, and GitHub, https://github.com/3DGenomes/TADkit

- `HiC-3DViewer` - HiC-3DViewer is an interactive web-based tool designed to provide an intuitive environment for investigators to facilitate the 3D exploratory analysis of Hi-C data. It based on Flask, it can be run directly or as a docker container. Bitbucket: https://bitbucket.org/nadhir/hic3dviewer/src/master/.
    - Mohamed Nadhir, Djekidel, Wang, Mengjie, Michael Q. Zhang, Juntao Gao. “HiC-3DViewer: a new tool to visualize Hi-C data in 3D space.” Quantitative Biology (2017) 5: 183. https://doi.org/10.1007/s40484-017-0091-8.

## De novo genome scaffolding

- Tools for de novo genome assembly from Hi-C reads: https://omictools.com/assembly-scaffolding-1-category

- `3D-DNA` Hi-C genome assembler and its application/validation. Methods are in the supplemental. https://github.com/theaidenlab/3D-DNA
    - Dudchenko, Olga, Sanjit S. Batra, Arina D. Omer, Sarah K. Nyquist, Marie Hoeger, Neva C. Durand, Muhammad S. Shamim, et al. “De Novo Assembly of the Aedes Aegypti Genome Using Hi-C Yields Chromosome-Length Scaffolds.” Science (New York, N.Y.) 356, no. 6333 (07 2017): 92–95. https://doi.org/10.1126/science.aal3327.

- `bin3C` - resolving metagenome-assembled genomes from Hi-C data. Metagenomic assembly using SPAdes (http://cab.spbu.ru/software/spades/). Tested using simulated (Sim3C and MetaART) and real-life data. Performance metrics: adjusted mutual information, weighted Bcubed. Contact matrix where bins are contigs. Infomap method for clustering the whole-contig graph. Compared with ProxiMeta (Phase Genomics). https://github.com/cerebis/bin3C
    - DeMaere, Matthew Z., and Aaron E. Darling. “Bin3C: Exploiting Hi-C Sequencing Data to Accurately Resolve Metagenome-Assembled Genomes.” Genome Biology 20, no. 1 (December 2019). https://doi.org/10.1186/s13059-019-1643-1.

- `dnaTri` - genome scaffolding via probabilistic modeling using two constrains of Hi-C data - distance-dependent decay and cis-trans ratio. Using known chromosome scaffolds and de novo assembly. Naive Bayes classifier to distinguish chromosome-specific vs. on different chromosomes contigs. Average linkage clustering to assemble contigs into 23 groups of chromosomes. Completed 65 previously unplaced contigs. Data, http://my5c.umassmed.edu/triangulation/, code https://github.com/NoamKaplan/dna-triangulation
    - Kaplan, Noam, and Job Dekker. “High-Throughput Genome Scaffolding from in Vivo DNA Interaction Frequency.” Nature Biotechnology 31, no. 12 (December 2013): 1143–47. https://doi.org/10.1038/nbt.2768.

- `GRAAL` - Genome (Re)Assembly Assessing Likelihood - genome assembly from Hi-C data. Gaps in genome assembly that can be filled by scaffolding. Superior than Lachesis and dnaTri, which are sensitive to duplications, clustering they use to initially arrange the scaffolds, parameters, unknown reliability. A Bayesian approach, prior assumptions are that cis-contact probabilities follow a power-law decay and that counts in the interaction matrix are Poisson. Multiple genomic structures tested using MCMC (Multiple-Try Metropolis algorithm) to maximize the likelihood of data given a genomic structure. https://github.com/koszullab/GRAAL and the next version instaGRAAL that uses https://github.com/koszullab/instaGRAAL
    - Marie-Nelly, Hervé, Martial Marbouty, Axel Cournac, Jean-François Flot, Gianni Liti, Dante Poggi Parodi, Sylvie Syan, et al. “High-Quality Genome (Re)Assembly Using Chromosomal Contact Data.” Nature Communications 5 (December 17, 2014): 5695. https://doi.org/10.1038/ncomms6695.

- `HiCAssembler` - Hi-C scaffolding tool combining assembly using Hi-C data with scaffolds from regular sequencing (short or long sequencing). Uses strategies from LACHESIS and 3D-DNA. Visual adjustment of scaffolding errors. Automatic and manual misassembly correction.  https://github.com/maxplanck-ie/HiCAssembler
    - Renschler, Gina, Gautier Richard, Claudia Isabelle Keller Valsecchi, Sarah Toscano, Laura Arrigoni, Fidel Ramirez, and Asifa Akhtar. “Hi-C Guided Assemblies Reveal Conserved Regulatory Topologies on X and Autosomes despite Extensive Genome Shuffling.” BioRxiv, March 18, 2019. https://doi.org/10.1101/580969.

- `Lachesis` - a three-step genome scaffolding tool: 1) graph clustering of scaffolds to chromosome groups, 2) ordering clustered scaffolds (minimum spanning tree, reassembling longest-to-shortest branches), 3) assigning orientation (exact position and the decay of interactions). Duplications and repeat regions may be incorrectly ordered/oriented. Tested on normal human, mouse, drosophila genomes, and on HeLa cancer genome. https://github.com/shendurelab/LACHESIS
    - Burton, Joshua N., Andrew Adey, Rupali P. Patwardhan, Ruolan Qiu, Jacob O. Kitzman, and Jay Shendure. “Chromosome-Scale Scaffolding of de Novo Genome Assemblies Based on Chromatin Interactions.” Nature Biotechnology 31, no. 12 (December 2013): 1119–25. https://doi.org/10.1038/nbt.2727.


## 3D reconstruction

- 3D genome reconstruction review. Intro into equilibrium/fractal globule models. Classification of reconstruction methods: distance-, contact-. and probability-based. [Table 1](https://biologicalproceduresonline.biomedcentral.com/articles/10.1186/s12575-019-0094-0#Tab1) summarizes many tools, methods, and references.
    - Oluwadare, Oluwatosin, Max Highsmith, and Jianlin Cheng. “An Overview of Methods for Reconstructing 3-D Chromosome and Genome Structures from Hi-C Data.” Biological Procedures Online 21, no. 1 (December 2019): 7. https://doi.org/10.1186/s12575-019-0094-0.

- `CSynth.org` - 3D genome browser visualization, GPU-based. Takes in text-based data in pair, BED, matrix, WIG, XYZ formats. 3D structure modeling, overlay of colored genomic annotations. http://csynth.org/
    - Todd, Stephen, Peter Todd, Simon J McGowan, James R Hughes, Yasutaka Kakui, Frederic Fol Leymarie, William Latham, and Stephen Taylor. “CSynth: A Dynamic Modelling and Visualisation Tool for 3D Chromatin Structure.” BioRxiv, January 3, 2019. https://doi.org/10.1101/499806.

- `Hierarchical3DGenome` - high-resolution (5kb) reconstruction of the 3D structure of the genome. Using LorDG (https://github.com/BDM-Lab/LorDG), first, assemble the 3D model at the level of TADs, then inside individual TADs. Gm12878 cell line, Arrowhead for TAD calling, KR and ICE normalization, benchmarking against miniMDS, five tests including comparison with FISH. https://github.com/BDM-Lab/Hierarchical3DGenome
    - Trieu, Tuan, Oluwatosin Oluwadare, and Jianlin Cheng. “Hierarchical Reconstruction of High-Resolution 3D Models of Large Chromosomes.” Scientific Reports 9, no. 1 (March 21, 2019): 4971. https://doi.org/10.1038/s41598-019-41369-w.

- `GenomeFlow` - a complete set of tools for Hi-C data alignment, normalization, 2D visualization, 3D genome modeling and visualization. ClusterTAD for TAD identification. LorDG and 3DMax for 3D genome reconstruction. https://github.com/jianlin-cheng/GenomeFlow
    - Trieu, Tuan, Oluwatosin Oluwadare, Julia Wopata, and Jianlin Cheng. “GenomeFlow: A Comprehensive Graphical Tool for Modeling and Analyzing 3D Genome Structure.” Bioinformatics (Oxford, England), September 12, 2018. https://doi.org/10.1093/bioinformatics/bty802.

- `ShRec3D` - shortest-path reconstruction in 3D. Genome reconstruction by translation a Hi-C matrix into a distance matrix, then multidimensional scaling. Uses binary contact maps. https://sites.google.com/site/julienmozziconacci/home/softwares
    - Lesne, Annick, Julien Riposo, Paul Roger, Axel Cournac, and Julien Mozziconacci. “3D Genome Reconstruction from Chromosomal Contacts.” Nature Methods 11, no. 11 (November 2014): 1141–43. https://doi.org/10.1038/nmeth.3104.


## Papers

A four-cutter enzyme yields a resolution of ∼256 bp and a six-cutter a resolution of ∼4,096 bp.

### Methodological Reviews

- Ay, Ferhat, and William S. Noble. “Analysis Methods for Studying the 3D Architecture of the Genome.” Genome Biology 16 (September 2, 2015): 183. https://doi.org/10.1186/s13059-015-0745-7. - Hi-C technology and methods review. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7#Tab1) - list of tools. Biases, normalization, matrix balancing. Extracting significant contacts, obs/exp ratio, parametric (powerlaw, neg binomial, double exponential), non-parametric (splines). 3D enrichment. References. TAD identification, directionality index. Outlook, importance of comparative analysis

- Chang, Pearl, Moloya Gohain, Ming-Ren Yen, and Pao-Yang Chen. “Computational Methods for Assessing Chromatin Hierarchy.” Computational and Structural Biotechnology Journal 16 (2018): 43–53. https://doi.org/10.1016/j.csbj.2018.02.003. - Review of higher-order (chromatin conformation capture) and primary order (DNAse, ATAC) technologies and analysis tools. Table 1 - technology summaries. Table 2 - tool summaries. Inter-chromosomal calls using Binarized contact maps. Visualization. Primary order technologies - details and peak calling.

- Forcato, Mattia, Chiara Nicoletti, Koustav Pal, Carmen Maria Livi, Francesco Ferrari, and Silvio Bicciato. “Comparison of Computational Methods for Hi-C Data Analysis.” Nature Methods, June 12, 2017. https://doi.org/10.1038/nmeth.4325. - Hi-C processing and TAD calling tools benchmarking, [Table 1](https://www.nature.com/articles/nmeth.4325/tables/1), simulated (Lun and Smyth method) and real data. Notes about pluses and minuses of each tool. TAD reproducibility is higher than chromatin interactions, increases with larger number of reads. Consistent enrichment of TAD boundaries in CTCF, irrespectively of TAD caller. Hi-C replication is poor, just a bit more than random. Supplementary table 2 - technical details about each program, Supplementary Note 1 - Hi-C preprocessing tools, Supplementary Note 2 - TAD callers. Supplementary note 3 - how to simulate Hi-C data. Supplementary note 6 - how to install tools. https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf

- Nicoletti, Chiara, Mattia Forcato, and Silvio Bicciato. “Computational Methods for Analyzing Genome-Wide Chromosome Conformation Capture Data.” Current Opinion in Biotechnology 54 (December 2018): 98–105. https://doi.org/10.1016/j.copbio.2018.01.023. - 3C-Hi-C tools review, Table 1 lists categorizes main tools, Figure 1 displays all steps in technology and analysis (alignment, resolution, normalization, including accounting for CNVs, A/B compartments, TAD detection, visualization). Concise description of all tools.

- Pal, Koustav, Mattia Forcato, and Francesco Ferrari. “Hi-C Analysis: From Data Generation to Integration.” Biophysical Reviews, December 20, 2018. https://doi.org/10.1007/s12551-018-0489-1. - Hi-C technology, data, 3D structures, analysis, and tools. Technology improvement and increasing resolution. FASTQ processing steps ("Hi-C data analysis: from FASTQ to interaction maps" section), pipelines, finding minimum resolution, normalization. Downstream analysis: A/B compartment detection, TAD callers, Hierarchical TADs, interaction callers. Data formats (pairix, sparse matrix format, cool, hic, butlr, hdf5, pgl). Hi-C visualization tools. Table 2 - summary and comparison of all tools.https://link.springer.com/article/10.1007%2Fs12551-018-0489-1#Tab2

- Yardımcı, Galip Gürkan, and William Stafford Noble. “Software Tools for Visualizing Hi-C Data.” Genome Biology 18, no. 1 (December 2017). https://doi.org/10.1186/s13059-017-1161-y. - Hi-C technology, data, and visualization review. Suggestion about graph representation.

- Waldispühl, Jérôme, Eric Zhang, Alexander Butyaev, Elena Nazarova, and Yan Cyr. “Storage, Visualization, and Navigation of 3D Genomics Data.” Methods, May 2018. https://doi.org/10.1016/j.ymeth.2018.05.008. - Review of tools for visualization of 3C-Hi-C data, challenges, analysis (Table 1). Data formats (hic, cool, BUTLR, ccmap). Database to quickly access 3D data. Details of each visualization tool in Section 4

### General Reviews

- Bouwman, Britta A. M., and Wouter de Laat. “Getting the Genome in Shape: The Formation of Loops, Domains and Compartments.” Genome Biology 16 (August 10, 2015): 154. https://doi.org/10.1186/s13059-015-0730-1. - TAD/loop formation review. Convergent CTCF, cohesin, mediator, different scenarios of loop formation. Stability and dynamics of TADs. Rich source of references.

- Chakraborty, Abhijit, and Ferhat Ay. “The Role of 3D Genome Organization in Disease: From Compartments to Single Nucleotides.” Seminars in Cell & Developmental Biology 90 (June 2019): 104–13. https://doi.org/10.1016/j.semcdb.2018.07.005. - 3D genome structure and disease. Evolution of technologies from FISH to variants of chromatin conformation capture. Hierarchical 3D organization, Table 1 summarizes each layer and its involvement in disease.. Rearrangement of TADs/loops in cancer and other diseases. Specific examples of the biological importance of TADs, loops as means of distal communication.

- Yu, Miao, and Bing Ren. “The Three-Dimensional Organization of Mammalian Genomes.” Annual Review of Cell and Developmental Biology 33 (06 2017): 265–89. https://doi.org/10.1146/annurev-cellbio-100616-060531. - 3D genome structure review. The role of gene promoters, enhancers, and insulators in regulating gene expression. Imaging-based tools, all flavors of chromatin conformation capture technologies. 3D features - chromosome territories, topologically associated domains (TADs), association of TAD boundaries with with replication domains, CTCF binding, transcriptional activity, housekeeping genes, genome reorganization during mitosis. Use of 3D data to annotate noncoding GWAS SNPs. 3D genome structure change in disease.

- Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation.” Molecular Systems Biology 11, no. 12 (December 23, 2015): 852–852. doi:10.15252/msb.20156492. http://msb.embopress.org/content/msb/11/12/852.full.pdf - 3D genome organization parts. Well-written and detailed. References. Technologies: FISH, 3C. 4C, 5C, Hi-C, GCC, TCC, ChIA-PET. Typical resolution - 40bp to 1Mb. LADs - conserved, but some are cell type-specific. Chromosome territories. Cell type-specific. inter-chromosomal interactions may be important to define cell-specific interactions. A/B compartments identified by PCA. Chromatin loops, marked by CTCF and Cohesin binding, sometimes, with Mediator. Transcription factories

- Dekker, Job, Marc A. Marti-Renom, and Leonid A. Mirny. “Exploring the Three-Dimensional Organization of Genomes: Interpreting Chromatin Interaction Data.” Nature Reviews. Genetics 14, no. 6 (June 2013): 390–403. https://doi.org/10.1038/nrg3454. https://www.nature.com/articles/nrg3454 - 3D genome review. Chromosomal territories, transcription factories. Details of each 3C technology. Exponential decay of interaction frequencies. Box 2: A/B compartments (several Mb), TAD definition, size (hundreds of kb). TADs are largely stable, A/B compartments are tissue-specific. Adjacent TADs are not necessarily of opposing signs, may jointly form A/B compartments. Genes co-expression, enhancer-promoters interactions are confined to TADs. 3D modeling.

- Witten, Daniela M., and William Stafford Noble. “On the Assessment of Statistical Significance of Three-Dimensional Colocalization of Sets of Genomic Elements.” Nucleic Acids Research 40, no. 9 (May 2012): 3849–55. https://doi.org/10.1093/nar/gks012.

### Normalization

- Yaffe, Eitan, and Amos Tanay. “Probabilistic Modeling of Hi-C Contact Maps Eliminates Systematic Biases to Characterize Global Chromosomal Architecture.” Nature Genetics 43, no. 11 (November 2011): 1059–65. https://doi.org/10.1038/ng.947. - Sources of biases: 1) non-specific ligation (large distance between pairs); 2) length of each ligated fragments; 3) CG content and nucleotide composition; 4) Mappability. Normalization. Enrichment of long-range interactions in active promoters. General aggregation of active chromosomal domains. Chromosomal territories, high-activity and two low-activity genomic clusters

### TAD detection

- Rocha, Pedro P., Ramya Raviram, Richard Bonneau, and Jane A. Skok. “Breaking TADs: Insights into Hierarchical Genome Organization.” Epigenomics 7, no. 4 (2015): 523–26. https://doi.org/10.2217/epi.15.25. - Textbook overview of TADs in 3 pages with key references. 3D organization discovery using FISH, 3C, Hi-C. Discovery of A/B compartments (euchromatin, heterochromatin), TADs as regulatory units conserved even across syntenic regions in different organisms. TADs coordinate gene expression. TAD boundaries are not created equal. Examples of changes of TAD boundaries (Hox gene cluster, ES differentiation). Hierarchy of TADs.

- Crane, Emily, Qian Bian, Rachel Patton McCord, Bryan R. Lajoie, Bayly S. Wheeler, Edward J. Ralston, Satoru Uzawa, Job Dekker, and Barbara J. Meyer. “Condensin-Driven Remodelling of X Chromosome Topology during Dosage Compensation.” Nature 523, no. 7559 (July 9, 2015): 240–44. https://doi.org/10.1038/nature14450. - InsulationScore, https://github.com/dekkerlab/crane-nature-2015 - Insulation score to define TADs - sliding square along the diagonal, aggregating signal within it. This aggregated score is normalized, and binned into TADs, boundaries. See Methods and implementation at https://github.com/dekkerlab/crane-nature-2015. ICE normalized data. OK to analyze data at two different resolutions

- Dali, Rola, and Mathieu Blanchette. “A Critical Assessment of Topologically Associating Domain Prediction Tools.” Nucleic Acids Research 45, no. 6 (April 7, 2017): 2994–3005. doi:10.1093/nar/gkx145. - TAD definition, tools. Meta-TADs, hierarchy, overlapping TADs. HiCPlotter for visualization. Manual annotation as a gold standard. Sequencing depth and resolution affects things. Code, manual annotations

- Forcato, Mattia, Chiara Nicoletti, Koustav Pal, Carmen Maria Livi, Francesco Ferrari, and Silvio Bicciato. “Comparison of Computational Methods for Hi-C Data Analysis.” Nature Methods, June 12, 2017. https://doi.org/10.1038/nmeth.4325. - Hi-C processing and TAD calling tools benchmarking, Table 1, simulated (Lun and Smyth method) and real data. Notes about pluses and minuses of each tool. TAD reproducibility is higher than chromatin interactions, increases with larger number of reads. Consistent enrichment of TAD boundaries in CTCF, irrespectively of TAD caller. Hi-C replication is poor, just a bit more than random. Supplementary table 2 - technical details about each program, Supplementary Note 1 - Hi-C preprocessing tools, Supplementary Note 2 - TAD callers. Supplementary note 3 - how to simulate Hi-C data. Supplementary note 6 - how to install tools. https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf. - Code: https://bitbucket.org/mforcato/hictoolscompare.git

- Olivares-Chauvet, Pedro, Zohar Mukamel, Aviezer Lifshitz, Omer Schwartzman, Noa Oded Elkayam, Yaniv Lubling, Gintaras Deikus, Robert P. Sebra, and Amos Tanay. “Capturing Pairwise and Multi-Way Chromosomal Conformations Using Chromosomal Walks.” Nature 540, no. 7632 (November 30, 2016): 296–300. https://doi.org/10.1038/nature20158. - TADs organize chromosomal territories. Active and inactive TAD properties. Methods: Good mathematical description of insulation score calculations. Filter TADs smaller than 250kb. Inter-chromosomal contacts are rare, ~7-10%. Concatemers (more than two contacts) are unlikely.

- Zufferey, Marie, Daniele Tavernari, Elisa Oricchio, and Giovanni Ciriello. “Comparison of Computational Methods for the Identification of Topologically Associating Domains.” Genome Biology 19, no. 1 (10 2018): 217. https://doi.org/10.1186/s13059-018-1596-9. - Comparison of 22 TAD callers across different conditions. Callers are classified as linear score-based, statistical model-based, clustering, graph theory. Table 1 and additional file 1 summarizes each caller. The effect of data resolution, normalization, hierarchy. Test on Rao 2014 data, chromosome 6. ICE or LGF (local genomic feature,) normalization. Measure of Concordance (MoC) to compare TADs. CTCF/cohesin as a measure of biological significance. TopDom, HiCseg, CaTCH, CHDF are the top performers. R scripts, including for calculation MoC, https://github.com/CSOgroup/TAD-benchmarking-scripts

- "Hierarchical Regulatory Domain Inference from Hi-C Data" - presentation by Bartek Wilczyński about TAD detection, existing algorithms, new SHERPA and OPPA methods. [Video](https://simons.berkeley.edu/talks/bartek-wilczynski-03-10-16), [PDF](https://simons.berkeley.edu/sites/default/files/docs/4588/2016-03-10-simons-institute-wilczynski.pdf), [Web site](http://regulomics.mimuw.edu.pl/wp/), [GitHub](https://github.com/regulomics/) - SHERPA and OPPA code there.

### TAD prediction

- Bednarz, Paweł, and Bartek Wilczyński. “Supervised Learning Method for Predicting Chromatin Boundary Associated Insulator Elements.” Journal of Bioinformatics and Computational Biology 12, no. 06 (December 2014): 1442006. doi:10.1142/S0219720014420062. http://www.worldscientific.com/doi/pdf/10.1142/S0219720014420062 - Predicting TAD boundaries using training data, and making new predictions. Bayesian network (BNFinder method), random forest vs. basic k-means clustering, ChromHMM, cdBEST. Using sequence k-mers and ChIP-seq data from modENCODE for prediction - CTCF ChIP-seq performs best. Used Boruta package for feature selection. Bayesian network performs best. To read on their BNFinder method

### TAD dynamics

- Zheng, H., and Xie, W. (2019). The role of 3D genome organization in development and cell differentiation. Nat. Rev. Mol. Cell Biol. https://www.nature.com/articles/s41580-019-0132-4 - 3D structure of the genome and its changes during gametogenesis, embryonic development, lineage commitment, differentiation. Changes in developmental disorders and diseases. Chromatin compartments and TADs. Chromatin changes during X chromosome inactivation. Promoter-enhancer interactions established during development are accompanied by gene expression changes. Polycomb-mediated interactions may repress developmental genes. References to many studies.

### Spectral clustering

- Y. X Rachel Wang, Purnamrita Sarkar, Oana Ursu, Anshul Kundaje and Peter J. Bickel, "Network modelling of topological domains using Hi-C data", https://arxiv.org/abs/1707.09587.  - TAD analysis using graph theoretical (network-based) methods. Treats TADs as a "community" within the network. Shows that naive spectral clustering is generally ineffective, leaving gaps in the data. 

- Liu, Sijia, Pin-Yu Chen, Alfred Hero, and Indika Rajapakse. “Dynamic Network Analysis of the 4D Nucleome.” BioRxiv, January 1, 2018. https://doi.org/10.1101/268318. - Temporal Hi-C data analysis using graph theory. Integrated with RNA-seq data. Network-based approaches such as von Neumann graph entropy, network centrality, and multilayer network theory are applied to reveal universal patterns of the dynamic genome. Toeplitz normalization. Graph Laplasian matrix. Detailed statistics.

- Norton, Heidi K, Harvey Huang, Daniel J Emerson, Jesi Kim, Shi Gu, Danielle S Bassett, and Jennifer E Phillips-Cremins. “Detecting Hierarchical 3-D Genome Domain Reconfiguration with Network Modularity,” November 22, 2016. https://doi.org/10.1101/089011. - Graph theory for TAD identification. Louvain-like local greedy algorithm to maximize network modularity. Vary resolution parameter, hierarchical TAD identification. Hierarchical spatial variance minimization method. ROC analysis to quantify performance. Adjusted RAND score to quantify TAD overlap.

- Chen, Jie, Alfred O. Hero, and Indika Rajapakse. “Spectral Identification of Topological Domains.” Bioinformatics (Oxford, England) 32, no. 14 (15 2016): 2151–58. https://doi.org/10.1093/bioinformatics/btw221. - Spectral algorithm to define TADs. Laplacian graph segmentation using Fiedler vector iteratively. Toeplitz normalization to remove distance effect. Spectral TADs do not overlap with Dixon's, but better overlap with CTCF.

- Fotuhi Siahpirani, Alireza, Ferhat Ay, and Sushmita Roy. “A Multi-Task Graph-Clustering Approach for Chromosome Conformation Capture Data Sets Identifies Conserved Modules of Chromosomal Interactions.” Genome Biology 17, no. 1 (December 2016). https://doi.org/10.1186/s13059-016-0962-8 - Arboretum-Hi-C - a multitask spectral clustering method to identify differences in genomic architecture. Intro about the 3D genome organization, TAD differences and conservation. Assessment of different clustering approaches using different distance measures, as well as raw contacts. Judging clustering quality by enrichment in regulatory genomic signals (Histone marks, LADs, early vs. late replication timing, TFs like POLII, TAF, TBP, CTCF, P300, CMYC, cohesin components, LADs, replication timing, SINE, LINE, LTR) and by numerical methods (Davies-Bouldin index, silhouette score, others). Although spectral clustering on contact counts performed best, spectral + Spearman correlation was chosen. Comparing cell types identifies biologically relevant differences as quantified by enrichment. Peak counts or average signal within regions were used for enrichment. Data https://zenodo.org/record/49767, and Arboretum-HiC https://bitbucket.org/roygroup/arboretum-hic


## URLs

- 4DN portal blog, https://www.4dnucleome.org/outreach.html

- `3d-genome-processing-tutorial` - A 3D genome data processing tutorial for ISMB/ECCB 2017. https://github.com/hms-dbmi/3d-genome-processing-tutorial

- `hic-data-analysis-bootcamp` - Workshop on measuring, analyzing, and visualizing the 3D genome with Hi-C data, https://github.com/hms-dbmi/hic-data-analysis-bootcamp
