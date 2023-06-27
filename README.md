# Hi-C data analysis tools and papers 

![MIT License](https://img.shields.io/badge/license-MIT-brightgreen.svg) 
[![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

Tools are added by publication date, newest on top. Unpublished tools are listed at the end of each section. See [Hi-C data notes](https://github.com/mdozmorov/HiC_data) and [single-cell Hi-C notes](https://github.com/mdozmorov/scHiC_notes) for more. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDnotes](https://github.com/mdozmorov/MDnotes) for other data science and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Pipelines](#pipelines)
  - [QC, quality control](#qc-quality-control)
  - [Capture-C](#capture-c)
    - [Capture-C peaks](#capture-c-peaks)
  - [HiChIP](#hichip)
  - [4C](#4c)
- [Resolution improvement](#resolution-improvement)
  - [Simulation](#simulation)
- [Normalization](#normalization)
  - [CNV-aware normalization](#cnv-aware-normalization)
- [Reproducibility](#reproducibility)
- [AB compartments](#ab-compartments)
- [Peak/Loop callers](#peak-loop-callers)
- [Differential analysis](#differential-analysis)
- [TAD callers](#tad-callers)
  - [TAD detection, benchmarking](#tad-detection-benchmarking)
  - [Architectural stripes](#architectural-stripes)
  - [Differential, timecourse TAD analysis](#differential-timecourse-tad-analysis)
- [Prediction of 3D features](#prediction-of-3d-features)
- [SNP-oriented](#snp-oriented)
- [CNV and Structural variant detection](#cnv-and-structural-variant-detection)
- [Visualization](#visualization)
- [De novo genome scaffolding](#de-novo-genome-scaffolding)
- [3D modeling](#3d-modeling)
- [Deconvolution](#deconvolution)
- [Haplotype phasing](#haplotype-phasing)
- [Papers](#papers)
  - [Methodological Reviews](#methodological-reviews)
  - [General Reviews](#general-reviews)
  - [Technology](#technology)
    - [Micro-C](#micro-c)
    - [Multi-way interactions](#multi-way-interactions)
    - [Imaging](#imaging)
  - [Normalization](#normalization-1)
  - [Spectral clustering](#spectral-clustering)
- [Courses](#courses)
- [Labs](#labs)
- [Misc](#misc)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Pipelines

- <a name="fan-c">[FAN-C](https://github.com/vaquerizaslab/fanc)</a> - Python pipeline for Hi-C processing. Input - raw FASTQ (aligned using BWA or Bowtie2, artifact filtering) or pre-aligned BAMs. KR or ICE normalization. Analysis and Visualization (contact distance decay, A/B compartment detection, TAD/loop detection, Average TAD/loop profiles, saddle plots, triangular heatmaps, comparison of two heatmaps). Automatic or modular. Compatible with .cool and .hic formats. [Tweet1](https://twitter.com/vaquerizas_lab/status/1225011187668209664?s=20), [Tweet2](https://twitter.com/vaquerizasjm/status/1339937903473025027?s=20). [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02215-9/tables/1) - detailed comparison of 13 Hi-C processing tools <details>
    <summary>Paper</summary>
    Kruse, Kai, Clemens B. Hug, and Juan M. Vaquerizas. "FAN-C: A Feature-Rich Framework for the Analysis and Visualisation of Chromosome Conformation Capture Data"  https://doi.org/10.1186/s13059-020-02215-9 Genome Biology 21, no. 1 (December 2020) 
</details>

- <a name="v3-galaxy-hicexplorer">[v3 of the Galaxy HiCExplorer](https://hicexplorer.usegalaxy.eu/)</a> - Includes full analysis of Hi-C, Capture-C, scHi-C. Workflow-like description of tools/tasks for each data type. <details>
  <summary>Paper</summary>
    Wolff, Joachim, Leily Rabbani, Ralf Gilsbach, Gautier Richard, Thomas Manke, Rolf Backofen, and Björn A Grüning. "Galaxy HiCExplorer 3: A Web Server for Reproducible Hi-C, Capture Hi-C and Single-Cell Hi-C Data Analysis, Quality Control and Visualization"  https://doi.org/10.1093/nar/gkaa220 Nucleic Acids Research, (July 2, 2020)
</details>

- <a name="schicexplorer">[scHiCExplorer](https://github.com/joachimwolff/scHiCExplorer)</a> - set of command-line tools specifically designed for scHi-C data. [scHiCExplorer's documentation](https://schicexplorer.readthedocs.io/en/latest/). <details>
  <summary>Paper</summary>
    Wolff, Joachim, Leily Rabbani, Ralf Gilsbach, Gautier Richard, Thomas Manke, Rolf Backofen, and Björn A Grüning. "Galaxy HiCExplorer 3: A Web Server for Reproducible Hi-C, Capture Hi-C and Single-Cell Hi-C Data Analysis, Quality Control and Visualization"  https://doi.org/10.1093/nar/gkaa220 Nucleic Acids Research, (July 2, 2020)
</details>

- [Cooltools](https://github.com/open2c/cooltools) formal paper - the suite of computational tools for modular high-level analysis of processed Hi-C data in [cooler format](https://github.com/open2c/cooler), by the [Open2C](https://open2c.github.io/) group. Supercede [hiclib](https://github.com/mirnylab/hiclib-legacy). Normalization (ICE, including trans-chromosomal), interaction frequency vs. distance decay curves (including within chromosomal arms), A/B compartment analysis (GC content, gene density for orientation, saddle plots), TAD/loop identification (insulation score, HiCCUPS, can work with up to 100bp resolution data, Micro-C, aggregate analyses, on- and off-diagonal pileups). Built using Python and command line interface. [Demo/documentation Jupyter notebools](https://cooltools.readthedocs.io/en/latest/notebooks/) demonstrating custom visualization and analysis, [GitHub](https://github.com/open2c/open2c_examples) version. [Code for manuscript figures](https://github.com/open2c/open2c_vignettes/tree/main/cooltools_manuscript). [Quaich](https://github.com/open2c/quaich) snakemake pipeline for Hi-C post-processing using [cooltools](https://github.com/open2c/cooltools), [chromosight](https://github.com/koszullab/chromosight, [mustache](https://github.com/ay-lab/mustache), [coolpuppy](https://github.com/open2c/coolpuppy). <details>
  <summary>Paper</summary>
  Open2C, Nezar Abdennur, Sameer Abraham, Geoffrey Fudenberg, Ilya M. Flyamer, Aleksandra A. Galitsyna, Anton Goloborodko, Maxim Imakaev, Betul A. Oksuz, and Sergey V. Venev. “Cooltools: Enabling High-Resolution Hi-C Analysis in Python.” Preprint. Bioinformatics, November 1, 2022. https://doi.org/10.1101/2022.10.31.514564.
</details>

- <a name="cooler">[Cooler: Scalable Storage for Hi-C Data and Other Genomically-Labeled Arrays](https://doi.org/10.1093/bioinformatics/btz540)</a>
    - [cooler](https://github.com/mirnylab/cooler)  - file format for storing Hi-C matrices, sparse, hierarchical, multi-resolution. `cooler` Python package for data loading, aggregation, merging, normalization (balancing), viewing, exporting data. Together with ["pairs" text-based format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md), and hic, cooler is accepted by the 4D Nucleome Consortium DAC. [Documentation](https://cooler.readthedocs.io/en/latest/)
    - [cooltools](https://github.com/mirnylab/cooltools) - tools to work with .cool files, [Documentation](https://cooltools.readthedocs.io/en/latest/)
    - [hiclib](https://github.com/mirnylab/hiclib-legacy) - Python tools to QC, map, normalize, filter and analyze Hi-C data
    - [hic2cool](https://github.com/4dn-dcic/hic2cool) - Lightweight converter between hic and cool contact matrices.
    - [pairtools](https://github.com/mirnylab/pairtools) - tools for low-level processing of mapped Hi-C paired reads. [Documentation](https://pairtools.readthedocs.io/en/latest/index.html).
  <details>
  <summary>Paper</summary>
    Abdennur, Nezar, and Leonid Mirny. "Cooler: Scalable Storage for Hi-C Data and Other Genomically-Labeled Arrays" https://doi.org/10.1093/bioinformatics/btz540 Bioinformatics, January 1, 2020
</details>

- [Pairtools](https://github.com/open2c/pairtools/) - chromatin conformation capture (Hi-C and more) modular CLI tools (Python implementation and ecosystem), from the aligned sam/bam files (bwa, minimap, etc.). Tools include _parse_ , _sort_, _dedup_. _parse_ considers various ligation types and properly reports contact pairs. Outputs a tab-separated [.pairs](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md) file (with header) with contact coordinates and additional information, can be converted into a binned contact matrix using cooler. _parse2_ handles multiple ligation events (multi-way contacts, PORE-C and MC-3C technologies), better handles newer sequencing data. _dedup_ handles imperfect matches, can consider additional columns. Additional preprocessing tools: _flip_, _header_, _select_, _sample_, _merge_. Quality control tools: _scaling_ and _stats, calculate chromosome-specific decay rate of interaction frequencies with distance, additional stats. Protocol-specific tools (_restrict_ for annotating pairs by restriction fragments, _phase_ for annotating haplotype-resolved Hi-C, _filterbycov_ for cleaning up single-cell Hi-C). Integration with [cooler](https://github.com/open2c/cooler) and [cooltools](https://github.com/open2c/cooltools) software from [Open2C](https://github.com/open2c), the backbone of the [distiller](https://github.com/open2c/distiller-nf) pipeline. Table 1 - comparison with Chromap, Juicer, HiC-Pro, HiCExplorer, Fan-C, TADbit. Second fastest to Chromap. pip, conda installable. [Documentation](https://pairtools.readthedocs.io/en/latest/). <details>
  <summary>Paper</summary>
  Open2C, Nezar Abdennur, Geoffrey Fudenberg, Ilya M. Flyamer, Aleksandra A. Galitsyna, Anton Goloborodko, Maxim Imakaev, and Sergey V. Venev. “Pairtools: From Sequencing Data to Chromosome Contacts.” Preprint. Bioinformatics, February 15, 2023. https://doi.org/10.1101/2023.02.13.528389.
</details>

- <a name="hicexplorer">[HiCExplorer](https://github.com/deeptools/HiCExplorer/)</a> - set of programs to process, normalize, analyze and visualize Hi-C data, Python, .cool format, conversion utilities. [Documentation](https://hicexplorer.readthedocs.io/). <details>
  <summary>Paper</summary>
  Ramírez, Fidel, Vivek Bhardwaj, Laura Arrigoni, Kin Chung Lam, Björn A. Grüning, José Villaveces, Bianca Habermann, Asifa Akhtar, and Thomas Manke. "High-Resolution TADs Reveal DNA Sequences Underlying Genome Organization in Flies"  https://doi.org/10.1038/s41467-017-02525-w Nature Communications 9, no. 1 (December 2018) 
</details> 

- <a name="galaxy-hicexplorer">[Galaxy HiCExplorer](https://hicexplorer.usegalaxy.eu/)</a> - a web server for Hi-C data preprocessing, QC, visualization. [Docker container](https://github.com/deeptools/docker-galaxy-hicexplorer). <details>
  <summary>Paper</summary>
    Wolff, Joachim, Vivek Bhardwaj, Stephan Nothjunge, Gautier Richard, Gina Renschler, Ralf Gilsbach, Thomas Manke, Rolf Backofen, Fidel Ramírez, and Björn A. Grüning. "Galaxy HiCExplorer: A Web Server for Reproducible Hi-C Data Analysis, Quality Control and Visualization"  https://doi.org/10.1093/nar/gky504 Nucleic Acids Research 46, no. W1 (July 2, 2018). 
</details>

- <a name="gitar">[GITAR](https://github.com/Zhong-Lab-UCSD/HiCtool)</a> - full Hi-C pre-processing, normalization, TAD detection, and visualization. Python scripts wrapping other tools. Table 1 summarizes the functionality of existing tools. [Documentation](https://www.genomegitar.org). <details>
    <summary>Paper</summary>
    Calandrelli, Riccardo, Qiuyang Wu, Jihong Guan, and Sheng Zhong. “GITAR: An Open Source Tool for Analysis and Visualization of Hi-C Data.” Genomics, Proteomics & Bioinformatics 16, no. 5 (2018): 365–72. https://doi.org/10.1016/j.gpb.2018.06.006 
</details>

- <a name="hic-bench">[HiC-bench](https://github.com/NYU-BFX/hic-bench)</a> - complete pipeline for Hi-C data analysis. <details>
    <summary>Paper</summary>
    Lazaris, Charalampos, Stephen Kelly, Panagiotis Ntziachristos, Iannis Aifantis, and Aristotelis Tsirigos. "HiC-Bench: Comprehensive and Reproducible Hi-C Data Analysis Designed for Parameter Exploration and Benchmarking"  https://doi.org/10.1186/s12864-016-3387-6 BMC Genomics 18, no. 1 (December 2017)
</details>

- <a name="tadbit">[TADbit](https://github.com/3DGenomes/tadbit)</a> - TADbit is a complete Python library to deal with all steps to analyze, model and explore 3C-based data. With TADbit, the user can map FASTQ files to obtain raw interaction binned matrices (Hi-C like matrices), normalize and correct interaction matrices, identify and compare the Topologically Associating Domains (TADs), build 3D models from the interaction matrices, and finally, extract structural properties from the models. TADbit is complemented by [TADkit](#tadkit) for visualizing 3D models. <details>
    <summary>Paper</summary>
    Serra, François, Davide Baù, Mike Goodstadt, David Castillo, Guillaume J. Filion, and Marc A. Marti-Renom. "Automatic Analysis and 3D-Modelling of Hi-C Data Using TADbit Reveals Structural Features of the Fly Chromatin Colors"  https://doi.org/10.1371/journal.pcbi.1005665 PLoS Computational Biology 13, no. 7 (July 2017) 
</details>

- <a name="juicer">[Juicer](https://github.com/theaidenlab/juicebox/wiki/Download)</a> - Java full pipeline to convert raw reads into Hi-C maps, visualized in Juicebox. Calls domains, loops, CTCF binding sites. `.hic` file format for storing multi-resolution Hi-C data. <details>
    <summary>Paper</summary>
    Durand, Neva C., Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. "Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments"  https://doi.org/10.1016/j.cels.2016.07.002 Cell Systems 3, no. 1 (July 2016)   

    Rao, Suhas S. P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. "A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping"  https://doi.org/10.1016/j.cell.2014.11.021 Cell 159, no. 7 (December 18, 2014)  - Juicer analysis example. TADs defined by frequent interactions. Enriched in CTCF and cohesin members. Five domain types. A1 and A2 enriched in genes. Chr 19 contains 6th pattern B6. Enrichment in different histone modification marks. TADs are preserved across cell types. Yet, differences between Gm12878 and IMR90 were detected. Boundaries detection by scanning image. Refs to the original paper.
</details>

- <a name="hic-pro">[HiC-Pro](https://github.com/nservant/HiC-Pro)</a> - Python and command line-based optimized and flexible pipeline for Hi-C data processing. [hicpro2juicebox](https://github.com/nservant/HiC-Pro/blob/master/bin/utils/hicpro2juicebox.sh) tool to generate Juicebox-compatible files (requires juicebox_clt.jar). [Documentation](http://nservant.github.io/HiC-Pro/) <details>
    <summary>Paper</summary>
    Servant, Nicolas, Nelle Varoquaux, Bryan R. Lajoie, Eric Viara, Chong-Jian Chen, Jean-Philippe Vert, Edith Heard, Job Dekker, and Emmanuel Barillot. "HiC-Pro: An Optimized and Flexible Pipeline for Hi-C Data Processing."  https://doi.org/10.1186/s13059-015-0831-x Genome Biology 16 (December 1, 2015) - HiC pipeline, references to other pipelines, comparison. From raw reads to normalized matrices. Normalization methods, fast and memory-efficient implementation of iterative correction normalization (ICE). Data format. Using genotyping information to phase contact maps.
</details>

- <a name="hicdat">[HiCdat](https://github.com/MWSchmid/HiCdat)</a> - Hi-C processing pipeline and downstream analysis/visualization. Analyses: normalization, correlation, visualization, comparison, distance decay, PCA, interaction enrichment test, epigenomic enrichment/depletion. Consists of GUI tool for data preprocessing and R package for data analysis. <details>
    <summary>Paper</summary>
    Schmid, Marc W., Stefan Grob, and Ueli Grossniklaus. "HiCdat: A Fast and Easy-to-Use Hi-C Data Analysis Tool"  https://doi.org/10.1186/s12859-015-0678-x BMC Bioinformatics 16 (September 3, 2015) 
</details>


- <a name="hicup">[HiCUP](http://www.bioinformatics.babraham.ac.uk/projects/hicup/)</a>  - Perl-based pipeline, alignment only, output BAM files. [Documentation](https://www.bioinformatics.babraham.ac.uk/projects/hicup/read_the_docs/html/index.html). <details>
    <summary>Paper</summary>
    Wingett, Steven, Philip Ewels, Mayra Furlan-Magaril, Takashi Nagano, Stefan Schoenfelder, Peter Fraser, and Simon Andrews. "HiCUP: Pipeline for Mapping and Processing Hi-C Data"  https://doi.org/10.12688/f1000research.7334.1 F1000Research 4 (2015) - HiCUP pipeline, alignment only, removes artifacts (religations, duplicate reads) creating BAM files. Details about Hi-C sequencing artifacts. Used in conjunction with other pipelines.
</details>

- <a name="hicpipeline">[HiC_Pipeline](https://github.com/XiaoTaoWang/HiC_pipeline)</a> - Python-based pipeline performing mapping, filtering, binning, and ICE-correcting Hi-C data, from raw reads (.sra, .fastq) to contact matrices. Additionally, converting to sparse format, performing QC. [Documentation](http://xiaotaowang.github.io/HiC_pipeline/).

- <a name="juicer-tools">[Juicer Tools](https://github.com/aidenlab/juicer/wiki/Juicer-Tools-Quick-Start)</a> - for creating/extracting data from .hic files, [Arrowhead](#arrowhead) for finding contact domains, [HiCCUPS](#hiccups) for loop detection and [HiCCUPS Diff](#hiccupsdiff) for finding differential loops, [MotifFinder](https://github.com/aidenlab/juicer/wiki/MotifFinder) for characterizing CTCF peaks at loop anchors, [Pearsons](https://github.com/aidenlab/juicer/wiki/Pearsons) for calculating the Pearson's correlation matrix of the Observed/Expected, [APA](https://github.com/aidenlab/juicer/wiki/APA) aggregated peak analysis, [Eigenvector](#eigenvector) for determining A/B chromatin states, [Compare Lists](https://github.com/aidenlab/juicer/wiki/Compare-Lists) for separating loop lists into common and condition-specific loops. All tools are described in the Extended Experimental Procedures of [Rao, Huntley et al. Cell 2014](https://www.cell.com/cms/10.1016/j.cell.2014.11.021/attachment/d3c6dcd8-c799-4f68-bbe4-201be54960b5/mmc1.pdf).

- [Juicer's 3D de novo assembly pipeline and resources](#theaidenlab3ddna).

- [ENCODE project Data Production and Processing Standard of the Hi-C Mapping Center, PDF](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf) - Computational standards of the Hi-C ENCODE mapping center including quality control measures and computational methods. 

- <a name="arima-mapping-pipeline">[Arima mapping pipeline](https://github.com/ArimaGenomics/mapping_pipeline)</a> - Mapping pipeline for data generated using Arima-HiC.

- <a name="cword">[cword](https://github.com/dekkerlab/cworld-dekker)</a> - Perl cworld module and collection of utility/analysis scripts for C data (3C, 4C, 5C, Hi-C).

- <a name="hicpipe">[HiCpipe](https://github.com/ChenFengling/HiCpipe)</a> - an efficient Hi-C data processing pipeline. It is based on [Juicer](#juicer) and [HiC-pro](#hic-pro), which combines the advantages of these two processing pipelines. HiCpipe is much faster than [Juicer](#juicer) and [HiC-Pro](#hic-pro) and can output multiple features of Hi-C maps.

- <a name="hicexperiment">[HiCExperiment](https://github.com/js2264/HiCExperiment)</a> - R package for handling three main Hi-C data formats ((m)cool, hic, HiC-Pro). Imports interaction pairs in the GInteractions objects with the intuitive metadata, like bin IDs, raw and normalized (balanced) interaction frequencies. The objects have the expected slots for features one can define from Hi-C data, like TADs, loops. Allows to query subsets (chunks) of Hi-C data.

- <a name="my5c">[my5C](http://my5c.umassmed.edu/)</a> - web-based tools, well-documented analysis and visualization of 5S data.

- <a name="nf-core-hic">[nf-core-hic](https://github.com/nservant/nf-core-hic)</a> - Analysis of Chromosome Conformation Capture data (Hi-C and more), [Nextflow](https://www.nextflow.io) pipeline. Also, [nf-core/hic](https://github.com/nf-core/hic). [Documentation](https://nf-co.re/hic/usage).

- <a name="distiller-nf">[distiller-nf](https://github.com/mirnylab/distiller-nf)</a> - Java modular Hi-C mapping pipeline for reproducible data analysis, [Nextflow](https://www.nextflow.io) pipeline. Alignment, filtering, aggregating Hi-C matrices.

- <a name="runhic">[runHiC](https://github.com/XiaoTaoWang/HiC_pipeline)</a> - aka `HiC_pipeline`, Hi-C data processing pipeline, from raw FASTQ files. Supports bwa-mem, chromap (default), and minimap2 aligners. Performs quality control (plots), filtering, binning (to mcool output), pileup. Parallelized. By [XiaoTao Wang](https://github.com/XiaoTaoWang).

- [4D Nucleome Hi-C Processing Pipeline](https://github.com/4dn-dcic/docker-4dn-hic) - set of scripts wrapped in a Docker image. Works with `.hic` and `.cool` files. [Overview](https://data.4dnucleome.org/resources/data-analysis/hi_c-processing-pipeline). 

- <a name="pairix">[Pairix](https://github.com/4dn-dcic/pairix)</a> - a tool to index and query files in [Pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md), a block-compressed text file format for storing paired genomic coordinates (header plus 7 columns: readID, chr1, pos1, chr2, pos2, strand1, strand2). Bgzipped sorted files (chr1, chr2, pos1, then pos2 sorting order) are indexed (less than a second for million lines) by pairix (similar in functionality, but incompatible with tabix). Command-line, R ([Rpairix](https://github.com/4dn-dcic/Rpairix)), Python implementations. Supplementary scripts like `bam2pairs`, `merged_nodups2pairs.pl`, `pairs_merger` are available. [pairsqc](https://github.com/4dn-dcic/pairsqc) - QC report generator for pairs files. Standard of the 4D Nucleome consortium, supported by [Juicer](https://github.com/aidenlab/juicer), [cooler](https://github.com/open2c/cooler), [pairtools](https://github.com/open2c/pairtools). <details>
    <summary>Paper</summary>
    Lee, Soohyun, Carl Vitzthum, Burak H. Alver, and Peter J. Park. "Pairs and Pairix: A File Format and a Tool for Efficient Storage and Retrieval for Hi-C Read Pairs"  https://doi.org/10.1101/2021.08.24.457552 Bioinformatics preprint, August 26, 2021.
</details>


### QC, quality control

- [qc3C](https://github.com/cerebis/qc3C) - Hi-C quality assessment method based on non-naturally occurring k-mers containing ligation artifacts (the proportion of "signal"). Details of various types for read-pairs, valid and invalid configurations. Tested on simulated and experimental data. Works of FASTQ or BAM files. Output - the breakdown of valid and invalid pairs (numbers, stacked barplots). Compatible with MultiQC. Conda, Docker, Singluarity installations. [Scripts for the paper](https://zenodo.org/record/4554522) <details>
    <summary>Paper</summary>
    DeMaere, Matthew Z., and Aaron E. Darling. "Qc3C: Reference-Free Quality Control for Hi-C Sequencing Data"  https://doi.org/10.1101/2021.02.24.432586 Preprint. Bioinformatics, February 25, 2021.
</details>

- <a name="hicnoisemeasurer">[HiCNoiseMeasurer](https://github.com/JRowleyLab/HiCNoiseMeasurer)</a> - a Python script to measure noise in .hic files using the auto-correlation function. 

- <a name="hicsampler">[HiCSampler](https://github.com/jrowleylab/hicsampler)</a> - a Python script for subsetting .hic files.

### Capture-C

- <a name="chicane">[CHiCANE](https://CRAN.R-project.org/package=chicane)</a> - an R-based data processing and interaction calling toolkit for the analysis and interpretation of Capture Hi-C data ([Arima](https://github.com/ArimaGenomics/CHiC), [Dovetail](https://github.com/dovetail-genomics/chicago)). Data preprocessing with [HiCUP](#hicup) (recommended), but BAMs from other pipelines are supported. Flexible regression modeling of the number of reads linking bait and target fragments, including distance. (Truncated) Negative binomial, Poisson distributions, zeros may be included, other covariates. Functionality to assess model fit. Similar tools - [GOTHiC](#gothic), [CHiCAGO](#chicago), [ChiCMaxima](#chicmaxima). Protocol explaining each step. <details>
    <summary>Paper</summary>
    Holgersen, Erle M. "Identifying High-Confidence Capture Hi-C Interactions Using CHiCANE" https://doi.org/10.1038/s41596-021-00498-1 NATURE PROTOCOLS (April 2021)
</details>

- [CaptureCompendium](http://userweb.molbiol.ox.ac.uk/public/telenius/CaptureCompendium/) - all-in-one toolkit for the design, analysis and presentation of 3C experiments, combines oligonucleotide design [Capsequm2](http://apps.molbiol.ox.ac.uk/CaptureC/cgi-bin/CapSequm.cgi), sequence mapping and extraction [CCseqBasic](https://github.com/Hughes-Genome-Group/CCseqBasicS), statistical data presentation and distribution [CaptureCompare](https://github.com/Hughes-Genome-Group/CaptureCompare) with [Peaky](#peaky) integration, [CaptureSee](https://capturesee.molbiol.ox.ac.uk/). Allows for multi-way interactions (Tri-C). [Overview of previous tools doing parts](https://userweb.molbiol.ox.ac.uk/public/telenius/CaptureCompendium/vignette/index.html). <details>
    <summary>Paper</summary>
    Telenius, Jelena M., Damien J. Downes, Martin Sergeant, A. Marieke Oudelaar, Simon McGowan, Jon Kerry, Lars L.P. Hanssen, et al. "CaptureCompendium: A Comprehensive Toolkit for 3C Analysis"  https://doi.org/10.1101/2020.02.17.952572 Preprint. Bioinformatics, February 18, 2020. 
</details>

- <a name="gopher">[GOPHER](https://github.com/TheJacksonLaboratory/Gopher)</a> - Java app probe design for Capture Hi-C. All, or selected, promoters, or around GWAS hits. [Documentation](https://gopher.readthedocs.io/en/latest/). <details>
    <summary>Paper</summary>
    Hansen, Peter, Salaheddine Ali, Hannah Blau, Daniel Danis, Jochen Hecht, Uwe Kornak, Darío G. Lupiáñez, Stefan Mundlos, Robin Steinhaus, and Peter N. Robinson. "GOPHER: Generator Of Probes for Capture Hi-C Experiments at High Resolution"  https://doi.org/10.1186/s12864-018-5376-4 BMC Genomics 20, no. 1 (December 2019). 
</details>

- [capC-MAP](https://github.com/cbrackley/capC-MAP) - Capture-C analysis pipeline. Python and C++, run through a configuration file. Outputs bedGraph. Compared with [HiC-Pro](#hic-pro), better detects PCR duplicates, identifies more interactions. Normalization tuned for Capture-C data. [Documentation](https://capc-map.readthedocs.io/) <details>
    <summary>Paper</summary>
    Buckle, Adam, Nick Gilbert, Davide Marenduzzo, and Chris A Brackley. "capC-MAP: software for analysis of Capture-C data"  https://doi.org/10.1093/bioinformatics/btz480 Bioinformatics, 15 November 2019
</details>



#### Capture-C peaks

- Benchmarking of Capture Hi-C analysis pipelines. HiCUP and mHiC preprocessing, the multimapping read rescue doesn't in mHiC doesn't improve data quality. GOTHIC, CHiCAGO, CHiCANE, CHiCMaxima tools comparison, reproducibility, the proportion of bait-bait interactions, overlap with open chromatin, histone marks. GOTHIC may be too permissive, CHiCANE is the strictest, CHiCAGO and CHiCMaxima overall provide good quality results. <details>
  <summary>Paper</summary>
  Aljogol, Dina, I. Richard Thompson, Cameron S. Osborne, and Borbala Mifsud. “Comparison of Capture Hi-C Analytical Pipelines.” Frontiers in Genetics 13 (January 28, 2022): 786501. https://doi.org/10.3389/fgene.2022.786501.
</details>

- <a name="chicago">[CHiCAGO](https://bioconductor.org/packages/Chicago/)</a> - protocol for Capture Hi-C analysis. Introduction into 3C-based technologies, as compared with Hi-C, Statistical model for background noise estimation, normalization, weighted p-value correction. Comparison with other tools ([HiCapTools](#hicaptools), [CHiCMaxima](#chicmaxima), [CHiCANE](#chicane)), downstream analysis with [Peaky](#peaky), [Chicdiff](#chicdiff). Preprocessing with [HiCUP](#hicup), input files (Table 1), how to create auxillary files and set parameters for different restriction enzymes (R, Python scripts), QC, visualization. [CHiCAGO R package](https://bioconductor.org/packages/Chicago/), [chicagoTools](https://bitbucket.org/chicagoTeam/chicago/src/master/), [PCHiCdata R package](https://bitbucket.org/chicagoTeam/chicago/src/master/PCHiCdata/inst/extdata/). <details>
    <summary>Paper</summary>
    Freire-Pritchett, Paula, Helen Ray-Jones, Monica Della Rosa, Chris Q. Eijsbouts, William R. Orchard, Steven W. Wingett, Chris Wallace, Jonathan Cairns, Mikhail Spivakov, and Valeriya Malysheva. "Detecting Chromosomal Interactions in Capture Hi-C Data with CHiCAGO and Companion Tools"  https://doi.org/10.1038/s41596-021-00567-5 Nature Protocols, August 9, 2021. 
    </details> <details>
    <summary>Paper</summary>
    Cairns, Jonathan, Paula Freire-Pritchett, Steven W. Wingett, Csilla Várnai, Andrew Dimond, Vincent Plagnol, Daniel Zerbino, et al. "CHiCAGO: Robust Detection of DNA Looping Interactions in Capture Hi-C Data."  https://doi.org/10.1186/s13059-016-0992-2 Genome Biology 17, no. 1 (2016): 127.
</details>

- <a name="peaky">[Peaky](https://github.com/cqgd/pky)</a> - Bayesian sparse variable selection approach. The model proposes that for any given bait, the expected CHi-C signal at each prey fragment is expressed as a sum of contributions from a set of fragments directly contacting that bait. [Documentation](https://cqgd.github.io/pky/articles/introduction.html). <details>
    <summary>Paper</summary>
    Eijsbouts, Christiaan Q, Oliver S Burren, Paul J Newcombe, and Chris Wallace. "Fine Mapping Chromatin Contacts in Capture Hi-C Data"  https://doi.org/10.1186/s12864-018-5314-5 BMC Genomics 20, no. 1 (December 2019). 
</details>

- <a name="chicmaxima">[ChiCMaxima](https://github.com/yousra291987/ChiCMaxima)</a> - a pipeline for detection and visualization of chromatin loops in Capture Hi-C data. Loess smoothing combined with a background model to detect significant interactions. Compare with [GOTHiC](#gothic) and [CHiCAGO](#chicago). <details>
    <summary>Paper</summary>
    Ben Zouari, Yousra, Anne M Molitor, Natalia Sikorska, Vera Pancaldi, and Tom Sexton. "ChiCMaxima: A Robust and Simple Pipeline for Detection and Visualization of Chromatin Looping in Capture Hi-C"  https://doi.org/10.1186/s13059-019-1706-3 Genome Biology, 22 May 2019
</details>

- <a name="hicaptools">[HiCapTools](https://github.com/sahlenlab/HiCapTools)</a> - A command-line software package that can design sequence capture probes for targeted chromosome capture applications and analyze sequencing output to detect proximities involving targeted fragments. Two probes are designed for each feature while avoiding repeat elements and non-unique regions. The data analysis suite processes alignment files to report genomic proximities for each feature at restriction fragment level and is isoform-aware for gene features. Statistical significance of contact frequencies is evaluated using an empirically derived background distribution. <details>
    <summary>Paper</summary>
    Anandashankar Anil, Rapolas Spalinskas, Örjan Åkerborg, Pelin Sahlén; "HiCapTools: a software suite for probe design and proximity detection for targeted chromosome conformation capture applications"  https://doi.org/10.1093/bioinformatics/btx625 Bioinformatics, Volume 34, Issue 4, 15 February 2018
</details>


### HiChIP

- <a name="hichip-peak">[HiChIP-Peak](https://github.com/ChenfuShi/HiChIP_peaks)</a> - Command-line HiChIP peak caller, focus on peaks at re-ligation sites. Peak filtering, then negative binomial model. Differential peak analysis similar to [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html). <details>
    <summary>Paper</summary>
    Shi, Chenfu, Magnus Rattray, and Gisela Orozco. "HiChIP-Peaks: A HiChIP Peak Calling Algorithm"  https://doi.org/10.1093/bioinformatics/btaa202 Bioinformatics, Volume 36, Issue 12, 15 June 2020
</details>

- <a name="maps">[MAPS](https://github.com/ijuric/MAPS)</a> - Model-based Analysis of PLAC-seq and HiChIP data. A zero-truncated Poisson regression framework to explicitly remove systematic biases, normalize, call long-distance interactions. Compared with [hichipper](https://github.com/aryeelab/hichipper), identifies more long-range, biologically relevant interactions. Works with [Arima's HiChIP data](https://arimagenomics.com/products/hichip/) <details>
    <summary>Paper</summary>
    Juric, Ivan, Miao Yu, Armen Abnousi, Ramya Raviram, Rongxin Fang, Yuan Zhao, Yanxiao Zhang, et al. "MAPS: Model-Based Analysis of Long-Range Chromatin Interactions from PLAC-Seq and HiChIP Experiments"  https://doi.org/10.1371/journal.pcbi.1006982 April 15, 2019, 24.
</details>

- <a name="cid">[CID](https://groups.csail.mit.edu/cgs/gem/cid/)</a> - Chromatin Interaction Discovery, call chromatin interactions from ChIA-PET. Outperforms [ChIA-PET2](https://github.com/GuipengLi/ChIA-PET2), [MANGO pipelines](https://github.com/dphansti/mango), call more peaks than [HICCUPS](#hiccups), [hichipper](https://github.com/aryeelab/hichipper). Java implementation <details>
    <summary>Paper</summary>
    Guo, Yuchun, Konstantin Krismer, Michael Closser, Hynek Wichterle, and David K Gifford. "High Resolution Discovery of Chromatin Interactions"  https://doi.org/10.1093/nar/gkz051 Nucleic Acids Research, February 14, 2019. 
</details>

- <a name="hichip-pipeline">[HiChIP pipeline](https://hichip.readthedocs.io/)</a> - by [Dovetail Genomics](https://dovetailgenomics.com/services/epigenetics/dovetail-hichip/), technology and analysis steps.

### 4C

- <a name="pipe4c">[pipe4C](https://github.com/deLaatLab/pipe4C)</a> - 4C-seq processing pipeline, R implementation. <details>
    <summary>Paper</summary>
    Krijger, Peter H.L., Geert Geeven, Valerio Bianchi, Catharina R.E. Hilvering, and Wouter de Laat. "4C-Seq from Beginning to End: A Detailed Protocol for Sample Preparation and Data Analysis"  https://doi.org/10.1016/j.ymeth.2019.07.014 Methods 170 (January 2020)
</details>

- <a name="peakC">[peakC](https://github.com/deWitLab/peakC)</a> - an R package for non-parametric peak calling in 4C/Capture-c/PCHiC data. <details>
    <summary>Paper</summary>
    Geeven, Geert, Hans Teunissen, Wouter de Laat, and Elzo de Wit. "PeakC: A Flexible, Non-Parametric Peak Calling Package for 4C and Capture-C Data"  https://doi.org/10.1093/nar/gky443 Nucleic Acids Research 46, no. 15 (September 6, 2018)
</details>

- <a name="4cseqpipe">[4Cseqpipe](https://github.com/changegene/4Cseqpipe)</a> processing pipeline and a genome-wide 4C primer database. <details>
    <summary>Paper</summary>
    Werken, Harmen J. G. van de, Gilad Landan, Sjoerd J. B. Holwerda, Michael Hoichman, Petra Klous, Ran Chachik, Erik Splinter, et al. "Robust 4C-Seq Data Analysis to Screen for Regulatory DNA Interactions"  https://doi.org/10.1038/nmeth.2173 Nature Methods 9, no. 10 (October 2012) - 4C technology paper. Two different 4bp cutters to increase resolution. Investigation of beta-globin locus, interchromosomal interactions.
</details>

## Resolution improvement

- <a name="vehicle">[VEHiCLE](https://github.com/Max-Highsmith/VEHiCLE)</a> - a variational autoencoder (feature extraction, dimensionality reduction) and Generative Adversarial Network (maps low-dimensional vectors to Hi-C maps) for Hi-C resolution enhancement. Uses a combination of four loss functions: adversarial loss, variational loss, mean square error, and insulation score loss (interesting!). Intro into VAEs, GANs, loss functions. Uses GM12878, IMR90, K562, HMEC data. Compared using five metrics (similarity, reproducibility) against [HiCPlus](#hicplus), [DeepHiC](#deephic), [HiCSR](#hicsr), outperforms all. Improves TAD identification, 3D structure modeling. Python implementation. <details>
    <summary>Paper</summary>
    Highsmith, Max, and Jianlin Cheng. "VEHiCLE: A Variationally Encoded Hi-C Loss Enhancement Algorithm for Improving and Generating Hi-C Data"  https://doi.org/10.1038/s41598-021-88115-9 Scientific Reports 11, no. 1 (December 2021)
</details>


- <a name="hicres">[HiCRes](https://github.com/ClaireMarchal/HiCRes)</a> - resolution estimation, based on the linear dependence of 20th percentile of coverage and the window size used to access coverage. Includes preseq for estimating and predicting library complexity, [bowtie2](https://github.com/BenLangmead/bowtie2) and [HiCUP](#hicup) for estimating Hi-C-specific QC metrics. Relatively insensitive to enzyme of choice. Implemented as Docker/Singularity images. Requires significant computational resources, like 5 hours on 40 CPU cluster. <details>
    <summary>Paper</summary>
    Marchal, Claire, Nivedita Singh, Ximena Corso-Díaz, and Anand Swaroop. "HiCRes: A Computational Method to Estimate and Predict the Resolution of HiC Libraries"  https://doi.org/10.1101/2020.09.22.307967 Preprint. Bioinformatics, September 22, 2020
</details>

- <a name="hicsr">[HiCSR](https://github.com/PSI-Lab/HiCSR)</a> - enhancement of Hi-C contact maps using a Generative Adversarial Network trained to optimize a custom loss function (weighted adversarial loss, pixel-wise L1 loss, and a feature reconstruction loss). An increase in resolution refers to recovering additional Hi-C contacts, "saturating" downsampled and noisy Hi-C matrices, not increasing the number of pixels. Representation learning with autoencoder with several convolutional layers and skip connections, then using it for the generator to create new matrices with discriminator telling them fake or real. Compared with [HiCPlus](#hicplus), [HiCNN](#hicnn), [hicGAN](#hicgan), [DeepHiC](#deephic). Reproducibility is better using four metrics. Python3 PyTorch implementation. <details>
    <summary>Paper</summary>
    Dimmick, Michael C., Leo J. Lee, and Brendan J. Frey. "HiCSR: A Hi-C Super-Resolution Framework for Producing Highly Realistic Contact Maps"  https://doi.org/10.1101/2020.02.24.961714 Preprint. Genomics, February 25, 2020. 
</details>

- <a name="deephic">[DeepHiC](https://github.com/omegahh/DeepHiC)</a> - a web-based generative adversarial network (GAN) for enhancing Hi-C data. Does not change the bin size, enhances the content of Hi-C data. Reconstructs the content from \~1% of the original data. Outperforms [Boost-HiC](#boost-hic), [HiCPlus](#hicplus), [HiCNN](#hicnn). [Documentation](http://sysomics.com/deephic/). <details>
    <summary>Paper</summary>
    - Hong, Hao, Shuai Jiang, Hao Li, Cheng Quan, Chenghui Zhao, Ruijiang Li, Wanying Li, et al. "DeepHiC: A Generative Adversarial Network for Enhancing Hi-C Data Resolution"  https://doi.org/10.1371/journal.pcbi.1007287 PLOS Computational Biology, February 21, 2020
</details>

- <a name="hifi">[HIFI](https://github.com/BlanchetteLab/HIFI)</a> - Command-line tool for Hi-C Interaction Frequency Inference for restriction fragment-resolution analysis of Hi-C data. Sparsity is resolved by using dependencies between neighboring restriction fragments, with Markov Random Fields performing the best. Better resolves TADs and sub-TADs, significant interactions. CTCF, RAD21, SMC3, ZNF143 are enriched around TAD boundaries. Matrices normalized for fragment-specific biases. <details>
    <summary>Paper</summary>
    - Cameron, Christopher JF, Josée Dostie, and Mathieu Blanchette. "Estimating DNA-DNA Interaction Frequency from Hi-C Data at Restriction-Fragment Resolution"  https://doi.org/10.1186/s13059-019-1913-y Genome Biology, 14 January 2020
</details>

- <a name="hicgan">[hicGAN](https://github.com/kimmo1019/hicGAN)</a> - improving resolution (saturation) of Hi-C data using Generative Adversarial Networks. Generator - five inner residual blocks to fight vanishing gradient (each block has two convolutional layers and batch normalization) and an outer skip connection. The discriminator has three convolutional blocks. Evaluation metrics: MSE, signal-to-noise ratio, structure similarity index, chromatin loop score. Compared against [HiCPlus](#hicplus). Python, Tensorflow implementation. <details>
    <summary>Paper</summary>
    Liu, Qiao, Hairong Lv, and Rui Jiang. "HicGAN Infers Super Resolution Hi-C Data with Generative Adversarial Networks"  https://doi.org/10.1093/bioinformatics/btz317 Bioinformatics 35, no. 14 (July 15, 2019)
</details>

- <a name="hicnn">[HiCNN](http://dna.cs.miami.edu/HiCNN/)</a> - a computational method for resolution enhancement. A modification of the [HiCPlus](#hicplus) approach, using very deep (54 layers, five types of layers) convolutional neural network. A Hi-C matrix of regular resolution is transformed into the high-resolution but very sparse matrix, HiCNN predicts the missing values. Pearson and MSE evaluation metrics, overlap of [Fit-Hi-C](#fit-hi-c)-detected significant interactions - perform similar or slightly better than [HiCPlus](#hicplus). PyTorch implementation. <details>
    <summary>Paper</summary>
    Liu, Tong, and Zheng Wang. "HiCNN: A Very Deep Convolutional Neural Network to Better Enhance the Resolution of Hi-C Data"  https://doi.org/10.1093/bioinformatics/btz251 Bioinformatics, April 9, 2019
</details>

- <a name="boost-hic">[Boost-HiC](https://github.com/LeopoldC/Boost-HiC)</a> - infer fine-resolution contact frequencies in Hi-C data, performs well even on 0.1% of the raw data. TAD boundaries remain. Better than [HiCPlus](#hicplus). It can be used for differential analysis (comparison) of two Hi-C maps. <details>
    <summary>Paper</summary>
    Carron, Leopold, Jean-baptiste Morlot, Vincent Matthys, Annick Lesne, and Julien Mozziconacci. "Boost-HiC: Computational Enhancement of Long-Range Contacts in Chromosomal Contact Maps"  https://doi.org/10.1101/471607 November 18, 2018. 
</details>

- <a name="mhi-c">[mHi-C](https://github.com/keleslab/mHiC)</a> - recovering alignment of multi-mapped reads in Hi-C data. Generative model to estimate probabilities for each bin-pair originating from a given origin. Reproducibility of contact matrices (stratum-adjusted correlation), reproducibility and number of significant interactions are improved. Novel interactions. Enrichment of TAD boundaries in LINE and SINE repetitive elements. Multi-mapping is not sensitive to trimming. Read filtering strategy (Figure 1, supplementary figures are very visual). <details>
    <summary>Paper</summary>
    Zheng, Ye, Ferhat Ay, and Sunduz Keles. "Generative Modeling of Multi-Mapping Reads with MHi-C Advances Analysis of High Throughput Genome-Wide Conformation Capture Studies"  https://doi.org/10.1101/301705 October 3, 2018.
</details>

- <a name="hicplus">[HiCPlus](https://github.com/zhangyan32/HiCPlus)</a> - increasing resolution of Hi-C data using convolutional neural network, mean squared error as a loss function. Basically, smoothing parts of Hi-C image, then binning into smaller parts. Performs better than bilinear/biqubic smoothing. <details>
    <summary>Paper</summary>
    Zhang, Yan, Lin An, Ming Hu, Jijun Tang, and Feng Yue. "HiCPlus: Resolution Enhancement of Hi-C Interaction Heatmap"  https://doi.org/10.1038/s41467-018-03113-2 March 1, 2017. 
</details>

### Simulation

- <a name="freehi-c-v.2.0">[FreeHi-C v.2.0](https://github.com/keleslab/FreeHiC)</a> - simulation of realistic Hi-C matrices with user- or data-driven spike-ins. Spike-ins are introduced on read-level and converted to interaction frequency level. Benchmark of [HiCcompare](#hiccompare), [multiHiCcompare](#multihiccompare), [diffHiC](#diffhic), and [Selfish](#selfish). Assessment of FDR, power, significance order, PRC and AUROC, genomic properties. GM12878 and A549 replicates of experimental Hi-C data. Three simulation settings with varying background distribution of interaction frequencies, spike-in proportions, sequencing depth. Figure 5 - summary of performances for all methods and comparison types. Subjective top performers: [multiHiCcompare](#multihiccompare), [HiCcompare](#hiccompare), [diffHiC](#diffhic), [Selfish](#selfish).<details>
    <summary>Paper</summary>
    Zheng, Ye, Peigen Zhou, and Sündüz Keleş. "FreeHi-C Spike-in Simulations for Benchmarking Differential Chromatin Interaction Detection"  https://doi.org/10.1016/j.ymeth.2020.07.001 Methods, July 2020
</details>

- <a name="freehi-c">[FreeHi-C](https://github.com/keleslab/FreeHiC)</a> - Hi-C data simulation based on properties of experimental Hi-C data. Preserves A/B compartments, TADs, the correlation between replicated ([HiCRep](https://github.com/TaoYang-dev/hicrep)), significant interactions, improves power to detect differential interactions. Robust to sequencing depth changes. Tested on replicates of GM12878, A549 human cancer cells, malaria *P.falciparum*. Compared with poorly performing [Sim3C](https://github.com/cerebis/sim3C). [Simulated data](https://zenodo.org/record/3345896). Python3 implementation. <details>
    <summary>Paper</summary>
    Zheng, Ye, and Sündüz Keleş. "FreeHi-C Simulates High-Fidelity Hi-C Data for Benchmarking and Data Augmentation"  https://doi.org/10.1038/s41592-019-0624-3 Nature Methods 17, no. 1 (January 2020)
</details>

## Normalization

- [HiConfidence](https://github.com/victorykobets/HiConfidence) - Python tool for eliminating biases from the Hi-C data by downweighting chromatin contacts from low-quality (low-coverage) Hi-C replicates. For each replicate, calculate differential matrix and then the confidence for each pixel as the inverse pixel-wise difference divided by their mean, raised to the power of a tunable parameter k (Figure 2A). Used for correction for replicates' confidence in calculating TAD boundaries, intra-TAD densities, improves replicate reproducibility (stratum-adjusted correlation coefficient). Compared with multiHiCcompare, aids in differential analysis, compartment and TAD detection. Applied to D. melanogaster S2 cells, [GSE200078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200078), processed with distiller, pairtools, cooltools. <details>
  <summary>Paper</summary>
  Kobets, Victoria A, Sergey V Ulianov, Aleksandra A Galitsyna, Semen A Doronin, Elena A Mikhaleva, Mikhail S Gelfand, Yuri Y Shevelyov, Sergey V Razin, and Ekaterina E Khrameeva. “HiConfidence: A Novel Approach Uncovering the Biological Signal in Hi-C Data Affected by Technical Biases.” Briefings in Bioinformatics, February 9, 2023, bbad044. https://doi.org/10.1093/bib/bbad044.
</details>

- <a name="hicorr">[HiCorr](https://github.com/JinLabBioinfo/HiCorr)</a> - a method for correcting known (mappability, CG content) and unknown (visibility) biases in Hi-C maps (multiplicative effects, Methods). Easy Hi-C protocol allowing for low-input (\~100K cells) Hi-C (in vivo HindIII digestion, in situ proximity ligation, DpnII digestion after lysis and reverse crosslink, Methods). HiCorr outputs ratio matrixes representing enrichment of Hi-C signal, hence loops can be easily extracted. Recovers 65% of [HICCUPS](#hiccups) loops and more. Chromatin loops are better marks of cell identity than compartments and outperform eQTLs in defining neurological GWAS target genes. Human iPSCs, neural progenitors (NPCs), neurons, fetal cerebellum, adult temporal cortex, data from other studies. <details>
    <summary>Paper</summary>
    Lu, Leina, Xiaoxiao Liu, Wei-Kai Huang, Paola Giusti-Rodríguez, Jian Cui, Shanshan Zhang, Wanying Xu, et al. "Robust Hi-C Maps of Enhancer-Promoter Interactions Reveal the Function of Non-Coding Genome in Neural Development and Diseases"  https://doi.org/10.1016/j.molcel.2020.06.007 Molecular Cell, June 2020
</details>

- <a name="normgam">[normGAM](http://dna.cs.miami.edu/normGAM/)</a> - an R package to normalize Genomic Architecture Mapping (GAM) data. New type of systematic bias, the fragment length bias. normGAM eliminates the fragment length bias resulting from random slicing, and biases related to window detection frequency, mappability, and GC content. Five normalization methods, including newly designed KR2 that handles negative values (others include the original GAM normalization algorithm normalized linkage disequilibrium NLD, VC, SCN, ICE). KR2 normalization produces better correlation with Hi-C data, all normalization methods improve concordance with FISH-detected distances.<details>
    <summary>Paper</summary>
    Liu, Tong, and Zheng Wang. "NormGAM: An R Package to Remove Systematic Biases in Genome Architecture Mapping Data"  https://doi.org/10.1186/s12864-019-6331-8 BMC Genomics, (December 2019)
</details>

- <a name="multihiccompare-norm">[multiHiCcompare](https://bioconductor.org/packages/multiHiCcompare/)</a> - R/Bioconductor package for joint normalization of multiple Hi-C datasets using cyclic loess regression through pairs of MD plots (minus-distance). Data-driven normalization accounting for the between-dataset biases. Per-distance [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)-based testing of significant interactions. <details>
    <summary>Paper</summary>
    Stansfield, John C, Kellen G Cresswell, and Mikhail G Dozmorov. "MultiHiCcompare: Joint Normalization and Comparative Analysis of Complex Hi-C Experiments"  https://doi.org/10.1093/bioinformatics/btz048 Bioinformatics, January 22, 2019. 
</details>

- <a name="binless">[Binless](https://github.com/3DGenomes/binless)</a> - a resolution-agnostic normalization method that adapts to the quality and quantity of available data, to detect significant interactions and differences. Negative binomial count regression framework, adapted for ICE normalization. Fused lasso to smooth neighboring signals. [TADbit](#tadbit) for data processing, details of read filtering. <details>
    <summary>Paper</summary>
    - Spill, Yannick G., David Castillo, Enrique Vidal, and Marc A. Marti-Renom. "Binless Normalization of Hi-C Data Provides Significant Interaction and Difference Detection Independent of Resolution"  https://doi.org/10.1038/s41467-019-09907-2 Nature Communications 10, no. 1 (26 2019)
</details>

- <a name="hiccompare_norm">[HiCcompare](https://bioconductor.org/packages/HiCcompare/)</a> - R/Bioconductor package for joint normalization of two Hi-C datasets using loess regression through an MD plot (minus-distance). Data-driven normalization accounting for the between-dataset biases. Per-distance permutation testing of significant interactions. <details>
    <summary>Paper</summary>
    Stansfield, John C., Kellen G. Cresswell, Vladimir I. Vladimirov, and Mikhail G. Dozmorov. "HiCcompare: An R-Package for Joint Normalization and Comparison of HI-C Datasets"  https://doi.org/10.1186/s12859-018-2288-x BMC Bioinformatics 19, no. 1 (December 2018). 
</details>

- <a name="hifive">[HiFive](https://www.taylorlab.org/software/hifive/)</a> - handling and normalization or pre-aligned Hi-C and 5C data. <details>
    <summary>Paper</summary>
    Sauria, Michael EG, Jennifer E. Phillips-Cremins, Victor G. Corces, and James Taylor. "HiFive: A Tool Suite for Easy and Efficient HiC and 5C Data Analysis"  https://doi.org/10.1186/s13059-015-0806-y Genome Biology 16, no. 1 (December 2015).  - HiFive - post-processing of aligned Hi-C and 5C data, three normalization approaches: "Binning" - model-based Yaffe & Tanay's method, "Express" - matrix-balancing approach, "Probability" - multiplicative probability model. Judging normalization quality by the correlation between matrices. 
</details>

- <a name="hicnorm">[HiCNorm](http://www.people.fas.harvard.edu/~junliu/HiCNorm/)</a> - removing known biases in Hi-C data (GC content, mappability, fragment length) via Poisson regression. <details>
    <summary>Paper</summary>
    Hu, Ming, Ke Deng, Siddarth Selvaraj, Zhaohui Qin, Bing Ren, and Jun S. Liu. "HiCNorm: Removing Biases in Hi-C Data via Poisson Regression"  https://doi.org/10.1093/bioinformatics/bts570 Bioinformatics (Oxford, England) 28, no. 23 (December 1, 2012) - Poisson normalization. Also tested negative binomial.
</details>

### CNV-aware normalization

- <a name="cancer-hic-norm">[Cancer-hic-norm](https://github.com/nservant/cancer-hic-norm)</a> - Hi-C data normalization considering CNVs. Extension of matrix-balancing algorithm to either retain the copy-number variation effect (LOIC) or remove them (CAIC). ICE itself can lead to misrepresentation of the contact probabilities between CNV regions. Estimating CNV directly from Hi-C data correcting for GC content, mappability, fragment length using Poisson regression. LOIC - the sum of contacts for a given genomic bin is proportional to CNV. CAIC - raw interaction counts are the product of a CNV bias matrix and the expected contact counts at a given genomic distance. [Data](http://members.cbio.mines-paristech.fr/~nvaroquaux/normalization/). LOIC and CAIC methods are implemented in the [iced](https://github.com/hiclib/iced) Python package. <details>
    <summary>Paper</summary>
    Servant, Nicolas, Nelle Varoquaux, Edith Heard, Emmanuel Barillot, and Jean-Philippe Vert. "Effective Normalization for Copy Number Variation in Hi-C Data"  https://doi.org/10.1186/s12859-018-2256-5 BMC Bioinformatics 19, no. 1 (September 6, 2018)
</details>

- <a name="oned">[OneD](https://github.com/qenvio/dryhic)</a> - CNV bias-correction method, addresses the problem of partial aneuploidy. Bin-centric counts are modeled using the negative binomial distribution, and its parameters are estimated using splines. A hidden Markov model is fit to infer the copy number for each bin. Each Hi-C matrix entry is corrected by dividing its value by the square root of the product of CNVs for the corresponding bins. Reproducibility score (eigenvector decomposition and comparison) to measure improvement in the similarity between replicated Hi-C data. <details>
    <summary>Paper</summary>
    Vidal, Enrique, François leDily, Javier Quilez, Ralph Stadhouders, Yasmina Cuartero, Thomas Graf, Marc A Marti-Renom, Miguel Beato, and Guillaume J Filion. "OneD: Increasing Reproducibility of Hi-C Samples with Abnormal Karyotypes"  https://doi.org/10.1093/nar/gky064 Nucleic Acids Research, January 31, 2018. 
</details>

- <a name="hicapp">[HiCapp](https://bitbucket.org/mthjwu/hicapp)</a> - Iterative correction-based caICB method. Method to adjust for the copy number variants in Hi-C data. Loess-like idea - we converted the problem of removing the biases across chromosomes to the problem of minimizing the differences across the count-distance curves of different chromosomes. Our method assumes equal representation of genomic locus pairs with similar genomic distances located on different chromosomes if there were no bias in the Hi-C maps. <details>
    <summary>Paper</summary>
    Wu, Hua-Jun, and Franziska Michor. "A Computational Strategy to Adjust for Copy Number in Tumor Hi-C Data"  https://doi.org/10.1093/bioinformatics/btw540 Bioinformatics (Oxford, England) 32, no. 24 (December 15, 2016)
</details>

## Reproducibility

- <a name="hprep">[HPRep](https://github.com/yunliUNC/HPRep)</a> - reproducibility measure for HiChIP and PLAC-seq data. HiCrep-inspired. Reorganize data into anchor-centric interaction bins, normalize (fragment length, GC content, mappability, log2 obs/exp) smooth, stratify by distance (concatenate bins with the same distance from anchors). Considering "AND" (bin-pairs), "XOR" (one anchor bin), "NOT" (no interactions, ignored) bin pairs. Distance metric - weighted Pearson correlation (pairs of columns) stratified by distance. Compared with [HiCRep](#hicrep), [HiC-Spector](#hic-spector), and naive Pearson on mouse H3K4me3 PLAC-seq data (brain and mESCs), human H3K37ac HiChIP data from GM12878 and K562, human H3K4me3 PLAC-seq brain data. HPRep shows higher similarity for replicates and more differentiation between cell lines, robust to downsampling. Nearly same results can be achieved analysing one chromosome (for speed). <details>
    <summary>Paper</summary>
    Rosen, Jonathan D, Yuchen Yang, Armen Abnousi, Jiawen Chen, Michael Song, Yin Shen, Ming Hu, and Yun Li. "HPRep: Quantifying Reproducibility in HiChIP and PLAC-Seq Datasets"  https://doi.org/10.3390/cimb43020082 Current issues in molecular biology, 17 September 2021
</details>

- [HiCrep.py](https://github.com/dejunlin/hicrep) - a fast Python implementation of stratum-adjusted correlation coefficient metric for measuring similarity between Hi-C datasets ([HiCrep](#hicrep) method, originally in R). Can be used for MDS. Evaluated on 90 datasets from 4D Nucleome. More than 20 times faster on a single CPU. Results are the same as R implementation.<details>
    <summary>Paper</summary>
    Lin, Dejun, Justin Sanders, and William Staﬀord Noble. "HiCRep.Py: Fast Comparison of Hi-C Contact Matrices in Python"  https://doi.org/10.1093/bioinformatics/btab097 Bioinformatics, February 12, 2021
</details>

- <a name="idr2d">[IDR2D](https://github.com/gifford-lab/idr2d)</a> - Irreproducible Discovery Rate that identifies replicable interactions in ChIP-PET, HiChIP, and Hi-C data. Includes the original 1D IDR version (https://github.com/nboley/idr). Resolves multiple pairwise interactions. <details>
    <summary>Paper</summary>
    Krismer, Konstantin, Yuchun Guo, and David K Gifford. "IDR2D Identifies Reproducible Genomic Interactions"  https://doi.org/10.1093/nar/gkaa030 Nucleic Acids Research, 06 April 2020
</details>

- [3DChromatin_ReplicateQC](https://github.com/kundajelab/3DChromatin_ReplicateQC) - Comparison of four Hi-C reproducibility assessment tools, [HiCRep](#hicrep), [GenomeDISCO](https://github.com/kundajelab/genomedisco), [HiC-Spector](#hic-spector), [QuASAR-Rep](https://github.com/kundajelab/3DChromatin_ReplicateQC). Tested the effects of noise, sparsity, resolution. Spearman doesn't work well. All tools performed similarly, worsening expectedly. [QuASAR](#quasar) has a QC tool measuring the level of noise. <details>
    <summary>Paper</summary>
    Yardimci, Galip, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. "Measuring the Reproducibility and Quality of Hi-C Data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1658-7 Genome Biology, March 19, 2019
</details>

- <a name="localtadism">[localtadsim](https://github.com/Kingsford-Group/localtadsim)</a> - Analysis of TAD similarity using a variation of information (VI) metric as a local distance measure. 23 human Hi-C datasets, [Hi-C Pro](#hic-pro) processed into 100kb matrices, [Armatus](#armatus) to call TADs. Defining structurally similar and variable regions. Comparison with previous studies of genomic similarity. Cancer-normal comparison - regions containing pan-cancer genes are structurally conserved in normal-normal pairs, not in cancer-cancer. <details>
    <summary>Paper</summary>
    - Sauerwald, Natalie, and Carl Kingsford. "Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types"  https://doi.org/10.1093/bioinformatics/bty265 Bioinformatics (Oxford, England) 34, no. 13 (July 1, 2018)
</details>

- <a name="quasar">[QuASAR](https://bxlab-hifive.readthedocs.io/en/latest/quasar_scoring.html)</a> - Hi-C quality and reproducibility measure using spatial consistency between local and regional signals. Finds the maximum useful resolution by comparing quality and replicate scores of replicates. Part of the [HiFive pipeline](#hifive). <details>
    <summary>Paper</summary>
    Sauria, Michael EG, and James Taylor. "QuASAR: Quality Assessment of Spatial Arrangement Reproducibility in Hi-C Data"  https://doi.org/10.1101/204438 BioRxiv, November 14, 2017.
</details>

- <a name="hicrep">[HiCRep](https://github.com/MonkeyLB/hicrep)</a> - Similarity assessment using generalized Cochran-Mantel-Haenzel statistics M2. Spearman/Pearson doesn't work. 2-step procedure: Smooth the matrix, then CMH statistics. Basically, splitting data by distance chunks, Pearson on each chunk, summarize. Simple and well-thought stats. Methods: Hi-C datasets with replicates, including 11 ENCODE datasets. [R package](https://github.com/MonkeyLB/hicrep), and [Python implementation](https://github.com/cmdoret/hicreppy) <details>
    <summary>Paper</summary>
    Yang, Tao, Feipeng Zhang, Galip Gurkan Yardimci, Ross C Hardison, William Stafford Noble, Feng Yue, and Qunhua Li. "HiCRep: Assessing the Reproducibility of Hi-C Data Using a Stratum-Adjusted Correlation Coefficient](https://genome.cshlp.org/content/27/11/1939.long Genome Research, August 30, 2017
</details>

- <a name="hic-spector">[HiC-Spector](https://github.com/gersteinlab/HiC-spector)</a> - reproducibility metric to quantify the similarity between contact maps using spectral decomposition. Decomposing Laplacian matrices and sum the Euclidean distance between eigenvectors. <details>
    <summary>Paper</summary>
    Yan, Koon-Kiu, Galip Gürkan Yardimci, Chengfei Yan, William S. Noble, and Mark Gerstein. "HiC-Spector: A Matrix Library for Spectral and Reproducibility Analysis of Hi-C Contact Maps"  https://doi.org/10.1093/bioinformatics/btx152 Bioinformatics (Oxford, England) 33, no. 14 (July 15, 2017)
</details>

## AB compartments

- <a name="pentad">[Pentad](https://github.com/magnitov/pentad)</a> - average A/B compartment analysis, within pre-defined range of genomic distances. A/B compartment areas from the observed-to-expected matrix rescaled using bilinear interpolation into squares of a predefined size, median averaged, plotted. Small (noisy) areas are filtered. Separate analysis of cis and trans interactions. Python implementation. Input: cooler Hi-C matrix and a compartment signal in the bedGraph format. <details>
  <summary>Paper</summary>
  Magnitov, Mikhail D., Azat K. Garaev, Alexander V. Tyakht, Sergey V. Ulianov, and Sergey V. Razin. “Pentad: A Tool for Distance-Dependent Analysis of Hi-C Interactions within and between Chromatin Compartments.” BMC Bioinformatics 23, no. 1 (December 2022): 116. https://doi.org/10.1186/s12859-022-04654-6.
</details>

- <a name="possum">[POSSUM](https://github.com/JRowleyLab/POSSUM)</a> - A/B compartment detection method in super-resolution Hi-C matrices. PCA of Sparse SUper Massive Matrices, Calculating eigenvectors for sparse matrices using power method (Figure 1, Methods). New GM12878 data at 500kb resolution (42 billion read pairs, 33 billion contacts). Genes can span compartments, but gene promoters almost exclusively (95%) are located in A compartments. Distinguishing loops formed by extrusion and non-extrusion mechanisms (SIP, [HiCCUPS](#hiccups), [Fit-Hi-C](#fit-hi-c) for detection), high resolution of Hi-C data is important. Applied to other datasets, organisms. A part of the Juicer pipeline [Eigenvector](#eigenvector), [C++ POSSUM code](https://github.com/JRowleyLab/POSSUM) on Jordan Rowley's lab GitHub. Other tools: [HiCSampler](https://github.com/jrowleylab/hicsampler), [HiCNoiseMeasurer](https://github.com/JRowleyLab/HiCNoiseMeasurer). [Tweet1](https://twitter.com/dphansti/status/1446126731723558915?s=20) by Doug Phanstiel, [Tweet2](https://twitter.com/MJordanRowley/status/1446127588032716810?s=20) by Jordan Rowley. <details>
    <summary>Paper</summary>
    - Gu, Huiya, Hannah Harris, Moshe Olshansky, Kiana Mohajeri, Yossi Eliaz, Sungjae Kim, Akshay Krishna, et al. "Fine-Mapping of Nuclear Compartments Using Ultra-Deep Hi-C Shows That Active Promoter and Enhancer Elements Localize in the Active A Compartment Even When Adjacent Sequences Do Not"  https://doi.org/10.1101/2021.10.03.462599 Preprint. Genomics, October 3, 2021.
</details>

- <a name="calder">[Calder](https://github.com/CSOgroup/CALDER)</a> - multi-scale compartment and sub-compartment detection, improvement over dichotomous AB compartment detection. Clustering contact similarities (Fisher's z-transformed correlations) into high intra and low inter-region similarities, followed by a divisive hierarchical clustering within each domain. The likelihood of nested sub-domains can be estimated using a mixture log-normal distribution. Detailed methods, complex. Eight subcompartments, 4 within the A and 4 within the B compartment, balanced set, in contrast to [SNIPER](#sniper). Expected associations with active/inactive genomic annotations. Nested compartments may be associated with TADs/loops. Analysis of domain repositioning across 114 cell lines. 40kb resolution. R package, named after [Alexander Calder](https://en.wikipedia.org/wiki/Alexander_Calder), an American sculptor. [Supplementary](https://www.nature.com/articles/s41467-021-22666-3#Sec23) Data 1 - IDs and links to Hi-C, ChIP-seq, and RNA-seq datasets; Data 2 - hg19 BED files of Complete domain hierarchies inferred by CALDER from 127 Hi-C contact maps; Data 7 - coordinates of Repositioned compartment domains between normal and cancer cell lines derived from breast, prostate, and pancreatic tissue samples. <details>
    <summary>Paper</summary>
    - Liu, Yuanlong, et al. "Systematic Inference and Comparison of Multi-Scale Chromatin Sub-Compartments Connects Spatial Organization to Cell Phenotypes"  https://doi.org/10.1038/s41467-021-22666-3 Nature Communication, 10 May 2021
</details>

- <a name="dchic">[dcHiC](https://github.com/ay-lab/dchic)</a> - differential A/B compartment analysis of Hi-C data. Uses Multiple Factor Analysis (MFA), and extension of PCA which combines Hi-C maps before performing generalized PCA. Analogous to weighted PCA in which every dataset is normalized for its biases (Methods). Multivariate distance measure to estimate statistical significance of compartment differences. Applied to mouse neuronal differentiation, mouse hematopoietic system, human cell Hi-C data. Gene enrichment analysis shows biologically relevant signal. Input - sparse matrix, hic, cool files. <details>
    <summary>Paper</summary>
    Wang, Jeffrey, Abhijit Chakraborty, and Ferhat Ay. "DcHiC: Differential Compartment Analysis of Hi-C Datasets"  https://doi.org/10.1101/2021.02.02.429297 BioRxiv, January 1, 2021
</details>

- <a name="sniper">[SNIPER](https://github.com/ma-compbio/SNIPER)</a> - 3D subcompartment (A1, A2, B1, B2, B3) identification from low-coverage Hi-C datasets. A neural network based on a denoising autoencoder (9 layers) and a multi-layer perceptron. Sigmoidal activation of inputs, ReLU, softmax on outputs. Dropout, binary cross-entropy. exp(-1/C) transformation of Hi-C matrices. Applied to Gm12878 and 8 additional cell types to compare subcompartment changes. Compared with Rao2014 annotations, outperforms Gaussian HMM and MEGABASE. <details>
    <summary>Paper</summary>
    Xiong, Kyle, and Jian Ma. "Revealing Hi-C Subcompartments by Imputing High-Resolution Inter-Chromosomal Chromatin Interactions"  https://doi.org/10.1038/s41467-019-12954-4 Nature Communications, 07 November 2019
</details>

- <a name="eigenvector">[Eigenvector](https://github.com/aidenlab/juicer/wiki/Eigenvector)</a> - Juicer's native tool. The eigenvector can be used to delineate compartments in Hi-C data at coarse resolution; the sign of the eigenvector typically indicates the compartment. The eigenvector is the first principal component of the Pearson's matrix.

## Peak/Loop callers

- Review of chromatin loop calling tools. Intro about loop formation (loop extrusion model), Hi-C data biases (GC content, mappability, fragment length, distance decay). Table 1 - loop calling tools for Hi-C, by year (model, language, application, etc.), Table 2 - loop calling for ChIA-PET, Table 3 - for HiChIP, Table 4 - for capture Hi-C. Table 5 - multiway interactions. Brief description of each method. <details>
  <summary>Paper</summary>
  Liu, Li, Kaiyuan Han, Huimin Sun, Lu Han, Dong Gao, Qilemuge Xi, Lirong Zhang, and Hao Lin. “A Comprehensive Review of Bioinformatics Tools for Chromatin Loop Calling.” Briefings in Bioinformatics 24, no. 2 (March 19, 2023): bbad072. https://doi.org/10.1093/bib/bbad072.
</details>

- [RefHiC](https://github.com/BlanchetteLab/RefHiC) - reference Hi-C data-guided TAD/loop detection (annotation). An attention-based deep learning frameworkthat determines which of the reference samples (4D Nucleome) are most relevant, and then makes a prediction based on the combined study sample and attention-weighted reference samples. Two components - a network combining the study sample and the reference panel and predicting loop points or left/right TAD boundary scores based on the local contact submatrix, and a task-specific component selecting one representative TAD/loop boundary (an encoder for dimensionality reduction, an attention module, a task-specific perceptron). Outperforms other tools (Mustache, Chromosight, HiCCUPS, Peakachu, RobusTAD and 13 TAD callers) across different cell types, species, and sequencing depths, using experimental ChIA-PET on CTCF, RAD21 and HiChIP on SMC1, H3K27ac. [Scripts to reproduce the paper](https://zenodo.org/record/7133194). Python. <details>
  <summary>Paper</summary>
  Zhang, Yanlin, and Mathieu Blanchette. “Reference Panel Guided Topological Structure Annotation of Hi-C Data.” Nature Communications 13, no. 1 (December 2, 2022): 7426. https://doi.org/10.1038/s41467-022-35231-3.
</details>

- <a name="lasca">[LASCA](https://github.com/ArtemLuzhin/LASCA_pipeline)</a> - loop/significant contact caller that uses Weibull distribution-based modeling to each diagonal. DBSCAN to cluster adjacent significant pixels. Works with Hi-C data from any species, tested on human, *C. Elegans*, *S. Cerevisiae*. Filters according Aggregate Peak Analysis patterns may be used to refine calls. Compared with [HiCCUPS](#hiccups), [MUSTACHE](#mustache), demonstrates good overlap. Also identifies non-CTCF-driven loops. Input - .cool files. [Python code](https://github.com/ArtemLuzhin/LASCA_pipeline). <details>
    <summary>Paper</summary>
    Luzhin, Artem V., Arkadiy K. Golov, Alexey A. Gavrilov, Artem K. Velichko, Sergey V. Ulianov, Sergey V. Razin, and Omar L. Kantidze. "LASCA: Loop and Significant Contact Annotation Pipeline"  https://doi.org/10.1038/s41598-021-85970-4 Scientific Reports, (December 2021)
</details>

- <a name="ziphi-c">[ZipHi-C](https://github.com/igosungithub/HMRFHiC)</a> - a Bayesian framework based on a Hidden Markov Random Field model to detect significant interactions and experimental biases in Hi-C data. Predecessors - [HMRFBayesHi-C](#hmrfbayeshic), [FastHiC](#fasthic). Borrows information from neighboring loci. Tested on simulated and experimental data, less false positives than [FastHiC](#fasthic), [Juicer](#juicer), [HiCExplorer](#hicexplorer). Detailed stats methods.<details>
    <summary>Paper</summary>
    Osuntoki, Itunu G., Andrew P. Harrison, Hongsheng Dai, Yanchun Bao, and Nicolae Radu Zabet. "ZipHiC: a novel Bayesian framework to identify enriched interactions and experimental biases in Hi-C data"  https://doi.org/10.1101/2021.10.19.463680 bioRxiv (October 20, 2021).
</details>

- <a name="cloops2">[cLoops2](https://github.com/YaqiangCao/cLoops2)</a> - improved pipeline for analysing Hi-TrAC/TrAC data. Peak/loop calling, differentially enriched loops, annotation, resolution estimation, similarity, aggregation/visualization. Improved blockDBSCAN clustering algorithm. Outperforms [MACS2](https://pypi.org/project/MACS2/), SICER, [HOMER](#homer), SEACR. <details>
    <summary>Paper</summary>
    Cao, Yaqiang, Shuai Liu, Gang Ren, Qingsong Tang, and Keji Zhao. "CLoops2: A Full-Stack Comprehensive Analytical Tool for Chromatin Interactions"  https://doi.org/10.1101/2021.07.20.453068 BioRxiv, 2021.
</details>

- <a name="neoloopfinder">[NeoLoopFinder](https://github.com/XiaoTaoWang/NeoLoopFinder)</a> - detecting chromatin interactions induced by all kinds of structural variants (SVs). Input - a Hi-C contact matrix and a list of SV breakpoints. Output - genome-wide CNV profile, CNV segments, local assembly around SVs (graph-based algorithm), corrected Hi-C matrix for newly assembled regions and normalized for CNV effect and allelic effect, chromatin loops in rearranged regions ([Peakachu](#peakachu)), enhancer-hijacking events (needs H3K27ac data). CNVs are detected by HMM-based segmentation module. Includes visualization module. Neo-loop detection in 50 cancer Hi-C datasets from cell lines and patient samples (17 cancer types). Cancer-specific neoloops, associated genes, epigenomic enrichments. Methods - DI + HMM. [Video, 20m](https://youtu.be/J61hFn5lB14). <details>
    <summary>Paper</summary>
    Wang, Xiaotao, Jie Xu, Baozhen Zhang, Ye Hou, Fan Song, Huijue Lyu, and Feng Yue. "Genome-Wide Detection of Enhancer-Hijacking Events from Chromatin Interaction Data in Rearranged Genomes"  https://doi.org/10.1038/s41592-021-01164-w Nature Methods, (June 2021)
</details>

- <a name="hic-act">[HiC-ACT](https://github.com/tmlagler/hicACT)</a> - improved chromatin loop detection considering spatial dependency (especially at high 5-10kb resolution). Aggregated Cauchy Test (ACT) based approach accounting for possible correlations between adjacent loci pairs from high-resolution Hi-C data. Combine a set of p-values, T statistics following Cauchy distribution under arbitrary dependence structure. Need the local smoothing bandwidth size. Post-processing of results from loop callers that assume independence among loci. Input - bin-pair identifiers and the corresponding p-values. Tested on GM12878 and mESC data. The improvement in power is most pronounced in low-depth (downsampled) data. Fast, implemented in R. [Documentation](https://yunliweb.its.unc.edu/hicACT/). <details>
    <summary>Paper</summary>
    Lagler, Taylor M., Armen Abnousi, Ming Hu, Yuchen Yang, and Yun Li. "HiC-ACT: Improved Detection of Chromatin Interactions from Hi-C Data via Aggregated Cauchy Test"  https://doi.org/10.1016/j.ajhg.2021.01.009 The American Journal of Human Genetics, (February 2021)
</details>

- <a name="chromosight">[Chromosight](https://github.com/koszullab/chromosight)</a> - python implemetnation of loop and pattern detection, computer vision-based (borders, FIREs, hairpins, and centromeres) in Hi-C maps. Takes in a single, whole-genome contact map, text-based bedGraph2d, and binary cool formats, ICE-normalizes. Sliding window, pattern detection using Pearson correlation with the template, then series of filters. Output - text-based. Outperforms [HiCexplorer](#hicexplorer), [HICCUPS](#hiccups), [HOMER](#homer), [cooltools](#cooler), in the order of decreasing F1. Tested on synthetic Hi-C data mimicking *S. cerevisiae* genome, [benchmark data at Zenodo](https://zenodo.org/record/3742095). <details>
    <summary>Paper</summary>
    Matthey-Doret, Cyril, Lyam Baudry, Axel Breuer, Rémi Montagne, Nadège Guiglielmoni, Vittore Scolari, Etienne Jean, et al. "Computer Vision for Pattern Detection in Chromosome Contact Maps"  https://doi.org/10.1038/s41467-020-19562-7 Nature Communications, (December 2020)
</details>

- <a name="peakachu">[Peakachu](https://github.com/tariks/peakachu)</a> - loop prediction from Hi-C data using random forest on loop-specific pixel intensities within 11x11 window. ChIA-PET and HiChIP provide positive training examples. H3K27ac HiChIP better predicts short-range interactions, CTCF ChIA-PET is better for longer interactions. 10Kb resolution data, Gm12878 and K562 cell line. Excels for short-range interactions. Detects more loops than [Fit-Hi-C](#fit-hi-c), [HiCCUPS](#hiccups), with good overlap. FDR estimated on auxin-treated vs. untreated HCT-116 cells, about 0.2%. Model trained using data from Hi-C performs well in other technologies, Micro-C, DNA SPRITE. Robust to sequencing depth. MCC to select best model. 3-fold cross-validation. Balanced training, same number of negative examples (with short and long distances between interacting loci). [Predicted loops for 56 cell/tissue types](http://promoter.bx.psu.edu/hi-c/publications.html). <details>
    <summary>Paper</summary>
    Salameh, Tarik J., Xiaotao Wang, Fan Song, Bo Zhang, Sage M. Wright, Chachrit Khunsriraksakul, Yijun Ruan, and Feng Yue. "A Supervised Learning Framework for Chromatin Loop Detection in Genome-Wide Contact Maps"  https://doi.org/10.1038/s41467-020-17239-9 Nature Communications, (December 2020)
</details>

- <a name="firecaller">[FIREcaller](https://yunliweb.its.unc.edu/FIREcaller/download.php)</a> - an R package to call frequently interacting regions from Hi-C data, as well as clustered super-FIREs. Normalization using HiCNormCis to regress out systematic biases. Converts normalized cis-interactions into Z-scores, calculates one-sided p-values and classifies bins as FIRE/nonFIRE. Also outputs continuous FIREscore (-ln(p-value)). FIREs are tissue-specific, can distinguish samples. Associated with H3K27ac and H3K4me3 signal. <details>
    <summary>Paper</summary>
    Crowley, Cheynna, Yuchen Yang, Yunjiang Qiu, Benxia Hu, Jakub Lipi, Hyejung Won, Bing Ren, Ming Hu, and Yun Li. "FIREcaller: Detecting Frequently Interacting Regions from Hi-C Data"  https://doi.org/10.1101/619288 October 26, 2020, 11.
</details>

- <a name="loopbit">[LOOPbit](https://github.com/3DGenomes/loopbit)</a> - loop detection guided by CTCF-CTCF topology classification. Interaction profiles of all CTCF-CTCF pairs are projected using self-organized feature maps (SOFM, [NeuPy](http://neupy.com/pages/home.html)), embedded using UMAP, clustered using HDBSCAN (10 clusters with strong-to-weak interaction patterns), and characterized by their epigenetic signatures (15 states aggregated into four). 10 clusters correspond to enrichment gradient from active to inactive chromatin states. These SOFM clusters serve as an input to a CNN. Trained CNN detects de novo chromatin loops from 9x9 Hi-C submatrices. Similar reproducibility as other loop callers, but epigenomically interpretable. Gm12878 data processed with [TADbit](#tadbit), [MetaWaffle](https://github.com/3DGenomes/metawaffle). Input - normalized Hi-C matrix in 3-column format. Python implementation, `scan` and `plot` functions. <details>
    <summary>Paper</summary>
    Galan, Silvia, François Serra, and Marc A. Marti-Renom. “Identification of Chromatin Loops from Hi-C Interaction Matrices by CTCF-CTCF Topology Classification.” Preprint. Bioinformatics, July 22, 2020. https://doi.org/10.1101/2020.07.21.214585.
</details>

- <a name="sip">[SIP](https://github.com/PouletAxel/SIP)</a> - loop caller using image analysis. Regional maxima-based, peaks called in a sliding window. Distance-normalized Hi-C matrices, image adjusted using Gaussian blur, contrast enhancement, White Top-Hat correction, identified peaks then filtered by peak enrichment, empirical FDR, loop decay. Comparison with [HiCCUPS](#hiccups) and [cLoops](#cloops) callers. Robust to noise, sequencing depth, much faster, good agreement, improved detection rate. SIPMeta - average metaplots of loops on bias-corrected images for better representation. Java implementation, works with .hic and .cool files. <details>
    <summary>Paper</summary>
    Rowley, M. Jordan, Axel Poulet, Michael H. Nichols, Brianna J. Bixler, Adrian L. Sanborn, Elizabeth A. Brouhard, Karen Hermetz, et al. "Analysis of Hi-C Data Using SIP Effectively Identifies Loops in Organisms from C. Elegans to Mammals"  https://doi.org/10.1101/gr.257832.119 Genome Research 30, no. 3 (March 2020)
</details>

- [HiCExplorer's hicDetectLoops](https://github.com/deeptools/HiCExplorer/) for loop detection - Review and critique of [HiCCUPS](#hiccups), [HOMER](#homer), [GOTHIC](#gothic), [cLoops](#cloops), [FastHiC](#fasthic). Distance-dependent of chromatin interactions with a continuous negative binomial distribution, detection of the interaction counts with p-values smaller than a threshold, then filtering. [Documentation](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html#hicdetectloops) <details>
    <summary>Paper</summary>
    Wolff, Joachim, Rolf Backofen, and Björn Grüning. "Loop Detection Using Hi-C Data with HiCExplorer"  https://doi.org/10.1101/2020.03.05.979096 Preprint. Bioinformatics, March 6, 2020
</details>

- <a name="mustache">[Mustache](https://github.com/ay-lab/mustache)</a> - loop detection from Hi-C and Micro-C maps. Scale-space theory, detection of blob-shaped objects in a multi-scale representation of contact maps, Gaussian kernels with increasing scales. Differences of adjacent Gaussians guide the search for local maxima. Series of filtering steps to minimize false positives. Corrected for multiple testing p-values of blobs. Applied to Gm12878 and K562 Hi-C data, and HFFc6 cell line Micro-C data, 5kb resolution. Compared with [HiCCUPS](#hiccups), SIP comparison added, detects similar and more loops flanked by convergent CTCF, RAD21, SMC3, loops confirmed by ChIA-PET and HiChIP data. Python3 tool, Conda/Docker wrapped, handles .hic/.cool files. [Tweet](https://twitter.com/ferhatay/status/1232716714346872832?s=20). <details>
    <summary>Paper</summary>
    Ardakany, Abbas Roayaei. "Mustache: Multi-Scale Detection of Chromatin Loops from Hi-C and Micro-C Maps Using Scale-Space Representation"  https://doi.org/10.1186/s13059-020-02167-0 Genome Biology, February 24, 2020
</details>

- <a name="fithic2">[FitHiC2](https://github.com/ay-lab/fithic)</a> - protocol to install/run FitHiC Python3 tool/scripts. Fit of non-increasing cubic splines to distance-interaction frequency decay to identify significant interactions in individual matrices. Accounts for biases derived from KR (ICE, or other) normalization (HiCKRy). Works with fixed-bin- or restriction cut site resolution data. Overview of FitHiC algorithm, accounting for biases. Flexible input options, from [HiC-Pro](#hic-pro), [Juicer](#juicer), and other tools, validPairs file format. Post-processing to prioritize highly significant interactions supported by the nearby loci, and filter noisy detections. HTML report, flexible BED-derived output format, conversion to formats for WashU epigenome and UCSC browsers. Installable using conda, pip, GitHub. Comparable methods - [HiCCUPS](#hiccups), [HOMER](#homer), [GOTHiC](#gothic), [HiC-DC](#hic-dc), a brief description of each. Tested on three datasets.Executable on [Code Ocean](https://codeocean.com/capsule/4528858/tree/v3), [Data](https://zenodo.org/record/3380589). <details>
    <summary>Paper</summary>
    Kaul, Arya, Sourya Bhattacharyya, and Ferhat Ay. "Identifying Statistically Significant Chromatin Contacts from Hi-C Data with FitHiC2"  https://doi.org/10.1038/s41596-019-0273-0 Nature Protocols, January 24, 2020.
</details>

- <a name="coolpup">[coolpup.py](https://github.com/Phlya/coolpuppy)</a> - Pile-up (aggregation, averaging) analysis of Hi-C data (.cool format) for visualizing and identifying chromatin loops from several sparse datasets, e.g., single-cell. Visualization using plotpup.py script. [Scripts for the paper](https://github.com/Phlya/coolpuppy_paper/tree/master/Nagano). <details>
    <summary>Paper</summary>
    Flyamer, Ilya M., Robert S. Illingworth, and Wendy A. Bickmore. "Coolpup.Py - a Versatile Tool to Perform Pile-up Analysis of Hi-C Data"  https://doi.org/10.1101/586537 BioRxiv, January 1, 2019
</details>

- <a name="cloops">[cLoops](https://github.com/YaqiangCao/cLoops)</a> - DBSCAN-based algorithm for the detection of chromatin loops in ChIA-PET, Hi-C, HiChIP, Trac-looping data. Local permutation-based estimation of statistical significance, several tests for enrichment over the background. Outperforms [diffHiC](#diffhic), [Fit-Hi-C](#fit-hi-c), [GOTHiC](#gothic), [HiCCUPS](#hiccups), [HOMER](#homer). <details>
    <summary>Paper</summary>
    Cao, Yaqiang, Xingwei Chen, Daosheng Ai, Zhaoxiong Chen, Guoyu Chen, Joseph McDermott, Yi Huang, and Jing-Dong J. Han. "Accurate Loop Calling for 3D Genomic Data with CLoops"  https://doi.org/10.1101/465849 November 8, 2018. 
</details>

- <a name="fithichip">[FitHiChIP](https://github.com/ay-lab/FitHiChIP)</a> - significant peak caller in HiChIP and PLAC-seq data. Accounts for assay-specific biases, as well as for the distance effect. 3D differential loops detection. Methods. <details>
    <summary>Paper</summary>
    Bhattacharyya, Sourya, Vivek Chandra, Pandurangan Vijayanand, and Ferhat Ay. "FitHiChIP: Identification of Significant Chromatin Contacts from HiChIP Data"  https://doi.org/10.1101/412833 September 10, 2018. 
</details>

- <a name="hic-dc">[HiC-DC](https://bitbucket.org/leslielab/hic-dc/src/master/)</a> - significant interaction detection using the zero-inflated negative binomial model and accounting for biases like GC content, mappability. Compared with Fit-Hi-C, more conservative. Robust to sequencing depth. Detects significant, biologically relevant interactions at all length scales, including sub-TADs. BWA-MEM alignment (Python script), then processing in R. <details>
    <summary>Paper</summary>
    Carty, Mark, Lee Zamparo, Merve Sahin, Alvaro González, Raphael Pelossof, Olivier Elemento, and Christina S. Leslie. "An Integrated Model for Detecting Significant Chromatin Interactions from High-Resolution Hi-C Data"  https://doi.org/10.1038/ncomms15454 Nature Communications 8, no. 1 (August 2017)
</details>

- <a name="gothic">[GoTHIC](https://www.bioconductor.org/packages/release/bioc/html/GOTHiC.html)</a> - R package for peak calling in individual HiC datasets, while accounting for noise. <details>
    <summary>Paper</summary>
    Mifsud, Borbala, Inigo Martincorena, Elodie Darbo, Robert Sugar, Stefan Schoenfelder, Peter Fraser, and Nicholas M. Luscombe. "GOTHiC, a Probabilistic Model to Resolve Complex Biases and to Identify Real Interactions in Hi-C Data"  https://doi.org/10.1371/journal.pone.0174744 Edited by Mark Isalan. PLOS ONE 12, no. 4 (April 5, 2017) - The GOTHiC (genome organization through HiC) algorithm uses a simple binomial distribution model to simultaneously remove coverage-associated biases in Hi-C data and detect significant interactions by assuming that the global background interaction frequency of two loci. Use of the Benjamini–Hochberg multiple-testing correction to control for the false discovery rate. 
</details>

- <a name="fasthic">[FastHiC](https://yunliweb.its.unc.edu/fasthic/)</a> - hidden Markov random field (HMRF)-based peak caller, fast and well-performing. <details>
    <summary>Paper</summary>
    Xu, Zheng, Guosheng Zhang, Cong Wu, Yun Li, and Ming Hu. "FastHiC: A Fast and Accurate Algorithm to Detect Long-Range Chromosomal Interactions from Hi-C Data"  https://doi.org/10.1093/bioinformatics/btw240 Bioinformatics (Oxford, England) 32, no. 17 (01 2016)
</details>

- <a name="hmrfbayeshic">[HMRFBayesHiC](https://yunliweb.its.unc.edu/HMRFBayesHiC/)</a> - a hidden Markov random field-based Bayesian peak caller to identify long-range chromatin interactions from Hi-C data. Borrowing information from neighboring loci. Previous peak calling methods, [Fit-Hi-C](#fit-hi-c). Interactions between enhancers and promoters as a benchmark. <details>
    <summary>Paper</summary>
    Xu, Zheng, Guosheng Zhang, Fulai Jin, Mengjie Chen, Terrence S. Furey, Patrick F. Sullivan, Zhaohui Qin, Ming Hu, and Yun Li. "A Hidden Markov Random Field-Based Bayesian Method for the Detection of Long-Range Chromosomal Interactions in Hi-C Data"  https://doi.org/10.1093/bioinformatics/btv650 Bioinformatics (Oxford, England) 32, no. 5 (01 2016)
</details>

- <a name="fit-hi-c">[Fit-Hi-C](https://noble.gs.washington.edu/proj/fit-hi-c/)</a> - Python tool for detection of significant chromatin interactions. <details>
    <summary>Paper</summary>
    Ay, Ferhat, Timothy L. Bailey, and William Stafford Noble. "Statistical Confidence Estimation for Hi-C Data Reveals Regulatory Chromatin Contacts"  https://doi.org/10.1101/gr.160374.113 Genome Research 24, no. 6 (June 2014) - Fit-Hi-C method, Splines to model distance dependence. Model mid-range interaction frequencies, decay with distance. Biases, normalization methods. Two-step splines - use all dots for the first fit, identify and remove outliers, second fit without outliers. Markers of boundaries - insulators, heterochromatin, pluripotent factors. CNVs are enriched in chromatin boundaries. Replication timing data how-to http://www.replicationdomain.com/. Validation Hi-C data. http://chromosome.sdsc.edu/mouse/hi-c/download.html
</details>

- <a name="hiccups">[HiCCUPS](https://github.com/aidenlab/juicer/wiki/HiCCUPS)</a> - chromatin loop detection, local maxima detection in Hi-C images. GPU and CPU implementations. Described in Section VI.a of the Extended Experimental Procedures of [Rao, Huntley et al. Cell 2014](https://www.cell.com/cms/10.1016/j.cell.2014.11.021/attachment/d3c6dcd8-c799-4f68-bbe4-201be54960b5/mmc1.pdf)

- <a name="hicpeaks">[HiCPeaks](https://github.com/XiaoTaoWang/HiCPeaks)</a> - Python CPU-based implementation for BH-FDR and [HICCUPS](#hiccups), two peak calling algorithms for Hi-C data, proposed by Rao et al. 2014. Text-to-[cooler](#cooler) Hi-C data converter, two scripts to call peaks, and one for visualization (creation of a .png file). [Pypi repo](https://pypi.org/project/hicpeaks/)

- <a name="homer">[HOMER](http://homer.ucsd.edu/homer/interactions/)</a> - Perl scripts for normalization, visualization, significant interaction detection, motif discovery. Does not correct for bias.

## Differential analysis

- Benchmarking of 11 methods for Hi-C map comparisons (Figure 1 outlines comparison goals and method classification). Basic methods (mean square error (MSE), Pearson/Spearman/Stratum-adjusted correlation, Structural similarity index), Map-informed methods (Eigenvector differences, directionality index, insulation differences, contact probability decay, triangle profiles), Feature-informed methods (cooltools TAD calling, HiCCUPS loop calling, overlap of functional genomics data). MSE and Spearman report different similarity properties (MSE is sensitive to intensity and Spearman is not). Directionality index performs best in identifying focal changes, like the loss of loops, while contact decay, eigenvector and insulation differences prioritize global changes. Basic methods are sensitive to technical variations (noise, resolution). Tested on experimental (ESC and HFF cells) and 22,500 in silico generated (unperturbed 1Mb sequences and with random 100bp structural variations, CTCF motif insertions/deletions, maps predicted with Akita) contact maps. The "Guidelines" section, Table 1, Supplementary Table 1 and Supplementary text - method description, strengths, weaknesses, and suggested applications of comparison methods. [GitHub](https://github.com/pollardlab/contact_map_scoring/) - data, Python code for all methods, Jupyter notebooks for all analyses. See also [HiC1Dmetrics](#hic1dmetrics) <details>
  <summary>Paper</summary>
  Gunsalus, Laura M, Evonne McArthur, Ketrin Gjoni, Shuzhen Kuang, John A Capra, and Katherine S Pollard. “Comparing Chromatin Contact Maps at Scale: Methods and Insights,” April 04, 2023, https://www.biorxiv.org/content/10.1101/2023.04.04.535480v1
</details>

- [StripeDiff](https://github.com/GuangyWang/stripeDiff) - differential stripe analysis. Intro about TADs/loops, stripes (observed at TAD boundaries, can be explained by loop extrusion), frequently interacting regions (FIREs, observed near the centers of TADs). Methods for stripe detection n single sample, fold change is typical. Definition of vertical/horizontal domains of a stripe window (by abrupt change around the window, Figure 1). Chi-square test for differences, AUROC for performance evaluation on simulated data (SimuStripe method) and manually annotated experimental data. Outperforms edgeR at various noise levels, stripes of different lengths, window sizes, sequencing depth, normalization strategies. Better detects stripes in single samples as compared to Zebra, StripeNN. Differential stripes are associated with CTCF binding changes, changes in chromatin state, gene expression, define cell lineage. [R code and data](https://zenodo.org/record/7026674) to reproduce the paper, [Python and R](https://github.com/GuangyWang/stripeDiff) implementation. <details>
  <summary>Paper</summary>
  Gupta, Krishan, Guangyu Wang, Shuo Zhang, Xinlei Gao, Rongbin Zheng, Yanchun Zhang, Qingshu Meng, Lili Zhang, Qi Cao, and Kaifu Chen. “StripeDiff: Model-Based Algorithm for Differential Analysis of Chromatin Stripe.” Science Advances 8, no. 49 (December 9, 2022): eabk2246. https://doi.org/10.1126/sciadv.abk2246.
</details>

- <a name="ne-mvnmf">[NE-MVNMF](https://github.com/Roy-lab/mvnmf)</a> - combines network enhancement (network denoising) and multitask non-negative matrix factorization for comparing multiple Hi-C matrices and idenrtifying dynamic (changing) regions. MVNMF - joint decomposition to find a common underlying structure in multiple matrices, clustering regions, comparing cluster assignments, identifying stretches of 5 or more regions switching clusters (dynamic regions). Applied to rat mammary epithelial cells, two time points (within and outside the window of susceptibility (WOS) to breast cancer, week 6 and 12). [Arima Genomics](https://arimagenomics.com), triplicates (6 samples total), merged into two, 10kb, [HiC-Pro](#hic-pro) processing, ICE normalization, differential point interactions as intersection of [Selfish](#selfish) and [Fit-Hi-C](#fit-hi-c). Integrated with gene expression (RNA-seq) and published breast cancer SNPs. [GEO GSE184285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184285). <details>
    <summary>Paper</summary>
    Baur, Brittany, Da-Inn Lee, Jill Haag, Deborah Chasman, Michael Gould, and Sushmita Roy. “Deciphering the Role of 3D Genome Organization in Breast Cancer Susceptibility.” Frontiers in Genetics 12 (January 11, 2022): 788318. https://doi.org/10.3389/fgene.2021.788318.
</details>

- <a name="hicdcplus">[HICDCPlus](https://bioconductor.org/packages/HiCDCPlus/)</a> - an R/Bioconductor package for Hi-C/Hi-ChIP interaction calling (directly from raw data, negative binomial regression accounting for genomic distance,GC content, mappability, restriction enzyme-based bin size) and differential analysis ([DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)). Includes TAD ([TopDom](#topdom)) and A/B compartment callers. Input - [HiC-Pro](#hic-pro) or [Juicer](#juicer). Output compatible with visualization in [Juicebox](#juicer) and [HiCExplorer](#hicexplorer). Compared with [diffHiC](#diffhic), [multiHiCcompare](#multihiccompare), [Selfish](#selfish) - provides better results. Normalization by ChIP-seq input may not be helpful. [BitBucket](https://bitbucket.org/leslielab/hicdcplus/src/master/). <details>
    <summary>Paper</summary>
    Sahin, Merve, Wilfred Wong, Yingqian Zhan, Kinsey Van Deynze, Richard Koche, and Christina S. Leslie. "HiC-DC+ enables systematic 3D interaction calls and differential analysis for Hi-C and HiChIP"  https://doi.org/10.1038/s41467-021-23749-x Nature communications, 07 June 2021
</details>

- <a name="bart3d">[BART3D](https://github.com/zanglab/bart3d)</a> - transcriptional regulators associated with differential chromatin interactions from Hi-C data. Input: [HiC-Pro](#hic-pro) matrices, .hic [Juicer](#juicer) files, .cool files. Output - ranked lists of transcriptional regulators. Distance-based normalization to average of individual matrices. Difference detection by a paired t-test of normalized interactions within 200kb (Figure 1A). Differential interactions are mapped to the union of DNAseI hypersensitive sites, then standard [BART](http://bartweb.org) algorithm. Python implementation. <details>
    <summary>Paper</summary>
    Wang, Zhenjia, Yifan Zhang, and Chongzhi Zang. "BART3D: Inferring Transcriptional Regulators Associated with Differential Chromatin InteracTions from Hi-C Data"  https://doi.org/10.1093/bioinformatics/btab173 Bioinformatics, 15 March 2021
</details>

- <a name="hic1dmetrics">[HiC1Dmetrics](https://github.com/wangjk321/HiC1Dmetrics)</a> - [Python3](https://pypi.org/project/h1d/), command-line and [web-based](http://hic1d.herokuapp.com/) tool for analyzing various types of 1D metrics for linerizing Hi-C data, identifying interesting 3D features (mostly TAD boundaries), comparing matrices. Review of 1D metrics, such as PC1, directionality index (DI), insulation score (IS), separation score (SS), contrast index (CI), distal-to-local ratio (DLR) and inter-chromosomal fraction of interactions (Figure 1, Methods). Newly developed intra-TAD score (IAS) and inter-TAD score (IES) (detection of stripes, asymmetric 5'/3' stripe TADs, loope, others), adjusted interaction frequency (IF) (detection of hubs, among others), directional relative frequency (DRF) for two-sample comparison, detects directional TADs (5'/3' dTADs) depicting an asymmetric event of inter-TAD interactions. New metrics can be compared between matrices, used for clustering, examples on experimental data. Novel observations, e.g., the upstream and downstream proximal TAD regions tend to exhibit opposite trends in gene expression. Inclusion of IS, PC1, and IF in ChromHMM identifies additional "active promoters on boundaries" state. Hi-C data preprocessing - Juicer, KR normalization. [Documentation](https://h1d.readthedocs.io). <details>
  <summary>Paper</summary>
  Wang, Jiankang, and Ryuichiro Nakato. “HiC1Dmetrics: Framework to Extract Various One-Dimensional Features from Chromosome Structure Data.” Briefings in Bioinformatics, December 1, 2021, bbab509. https://doi.org/10.1093/bib/bbab509.
</details>

- <a name="chess">[CHESS](https://github.com/vaquerizaslab/CHESS)</a> (Comparison of Hi-C Experiments using Structural Similarity) - comparative analysis of Hi-C matrices and automatic feature extraction (TADs, loops, stripes). Image analysis-based structural similarity index (SSIM, combines brightness, contrast, and structure differences, S = 1 - identical matrices, <1 - differences) to assign similarity score and an associated p-value (from empirical distribution of SSIMs, two types of bacground model, Methods) to pairs of genomic regions. Obs/exp transformation, differential matrix, denoise, smooth, binarize, feature extraction using close morphology filter, k-means clustering, classification. Works with low sequencing depth, high noise data. Outperforms [diffHiC](#diffhic), [HOMER](#homer), and [ACCOST](#accost). Applied to interspecies comparison of syntenic regions (Synteny portal), WT and Zelda-depleted Drosophila, C-cell lymphoma, Capture-C analysis. Input - [Juicer](#juicer), [Cooler](#cooler), of [FAN-C](#fan-c) format, plus .bedpe for regions to compare (chess pairs to generate). Python, scikit-image module. [Documentation](https://chess-hic.readthedocs.io/en/latest/index.html). <details>
    <summary>Paper</summary>
    Galan, Silvia, Nick Machnik, Kai Kruse, Noelia Díaz, Marc A. Marti-Renom, and Juan M. Vaquerizas. "CHESS Enables Quantitative Comparison of Chromatin Contact Data and Automatic Feature Extraction"  https://doi.org/10.1038/s41588-020-00712-y Nature Genetics, October 19, 2020
</details>

- <a name="diffgr">[DiffGR](https://github.com/wmalab/DiffGR)</a> - differentially interacting genomic regions. Stratum-adjusted correlation coefficient (SCC) (HiCrep-inspired) to measure similarity of local TAD regions. Focus on within-TAD interactions. Simulated data at various levels of sparsity, noise, [HiCseg](#hicseg) for TAD calling. 2D mean filter for smoothing, KR normalization. Permutation test to estimate the significance of SCC changes. FDR depends on the proportion of altered TADs. R implementation. <details>
    <summary>Paper</summary>
    Liu, Huiling, and Wenxiu Ma. "DiffGR: Detecting Differentially Interacting Genomic Regions from Hi-C Contact Maps"  https://doi.org/10.1101/2020.08.29.273698 bioRxiv, August 31, 2020
</details>

- <a name="serpentine">[Serpentine](https://github.com/koszullab/serpentine)</a> - differential analysis of two Hi-C maps using the 2D serpentine-binning method. Serpentine is a subset of connected pixels defined by thresholds in control and experimental contact maps. Serpentines are then compared using the Mean-Deviation plot. Help to alleviate the effect of sparsity. Uses [HiCcompare](#hiccompare) functionality. Normalization does not help. Python package, currently processes full 1500x1500 matrices. [Documentation](https://serpentine.readthedocs.io/en/latest/) <details>
    <summary>Paper</summary>
    Baudry, Lyam, Gaël A Millot, Agnes Thierry, Romain Koszul, and Vittore F Scolari. "Serpentine: A Flexible 2D Binning Method for Differential Hi-C Analysis"  https://doi.org/10.1093/bioinformatics/btaa249 Edited by Alfonso Valencia. Bioinformatics 36, no. 12 (June 1, 2020)
</details>

- <a name="accost">[ACCOST](https://bitbucket.org/noblelab/accost/src/master/)</a> (Altered Chromatin COnformation STatirstics) - distance-aware differential Hi-C analysis. Extends the statistical model of DEseq by using the size factors to model the genomic distance effect. Use of the MD plot. Compare with [diffHiC](#diffhic), [FIND](#find), and [HiCcompare](#hiccompare). Evaluated on human, mouse, plasmodium Hi-C data. <details>
    <summary>Paper</summary>
    - Cook, Kate B., Borislav H. Hristov, Karine G. Le Roch, Jean Philippe Vert, and William Stafford Noble. "Measuring significant changes in chromatin conformation with ACCOST"  https://doi.org/10.1093/nar/gkaa069 Nucleic acids research, (18 March 2020) 
</details>

- <a name="selfish">[Selfish](https://github.com/ucrbioinfo/Selfish)</a> - comparative analysis of replicate Hi-C experiments via a self-similarity measure - local similarity borrowed from image comparison. Check reproducibility, detect differential interactions. Boolean representation of contact matrices for reproducibility quantification. Deconvoluting local interactions with a Gaussian filter (putting a Gaussian bell around a pixel), then comparing derivatives between contact maps for each radius. Simulated (Zhou method) and real comparison with [FIND](#find) - better performance, especially on low fold-changes. Stronger enrichment of relevant epigenomic features. Matlab implementation. <details>
    <summary>Paper</summary>
    Roayaei Ardakany, Abbas, Ferhat Ay, and Stefano Lonardi. "Selfish: Discovery of Differential Chromatin Interactions via a Self-Similarity Measure"  https://doi.org/10.1093/bioinformatics/btz362 Bioinformatics, July 2019
</details>

- <a name="multihiccompare">[multiHiCcompare](https://bioconductor.org/packages/multiHiCcompare/)</a> - R/Bioconductor package for joint normalization of multiple Hi-C datasets using cyclic loess regression through pairs of MD plots (minus-distance). Data-driven normalization accounting for the between-dataset biases. Per-distance [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)-based testing of significant interactions. <details>
    <summary>Paper</summary>
    Stansfield, John C, Kellen G Cresswell, and Mikhail G Dozmorov. "MultiHiCcompare: Joint Normalization and Comparative Analysis of Complex Hi-C Experiments"  https://doi.org/10.1093/bioinformatics/btz048 Bioinformatics, January 22, 2019
</details>

- <a name="chicdiff">[Chicdiff](https://github.com/RegulatoryGenomicsGroup/chicdiff)</a> - differential interaction detection in Capture Hi-C data. Signal normalization based on the [CHiCAGO](#chicago) framework, differential testing using [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html). Accounting for distance effect by the Independent Hypothesis Testing (IHW) method to learn p-value weights based on the distance to maximize the number of rejected null hypotheses. <details>
    <summary>Paper</summary>
    Cairns, Jonathan, William R. Orchard, Valeriya Malysheva, and Mikhail Spivakov. "Chicdiff: A Computational Pipeline for Detecting Differential Chromosomal Interactions in Capture Hi-C Data"  https://doi.org/10.1101/526269 BioRxiv, January 1, 2019
</details>

- <a name="hiccompare">[HiCcompare](https://bioconductor.org/packages/HiCcompare/)</a> - R/Bioconductor package for joint normalization of two Hi-C datasets using loess regression through an MD plot (minus-distance). Data-driven normalization accounting for the between-dataset biases. Per-distance permutation testing of significant interactions. <details>
    <summary>Paper</summary>
    Stansfield, John C., Kellen G. Cresswell, Vladimir I. Vladimirov, and Mikhail G. Dozmorov. "HiCcompare: An R-Package for Joint Normalization and Comparison of HI-C Datasets"  https://doi.org/10.1186/s12859-018-2288-x BMC Bioinformatics 19, no. 1 (December 2018)
</details>

- <a name="find">[FIND](https://bitbucket.org/nadhir/find/src/master/)</a> - differential chromatin interaction detection comparing the local spatial dependency between interacting loci. Previous strategies - simple fold-change comparisons, binomial model ([HOMER](#homer)), count-based [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). FIND exploits a spatial Poisson process model to detect differential chromatin interactions that show a significant change in their interaction frequency and the interaction frequency of their adjacent bins. "Variogram" concept. For each point, compare densities between conditions using Fisher's test. Explored various multiple correction testing methods, used r^th ordered p-values (rOP) method. Benchmarking against [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) in simulated settings - FIND outperforms at shorter distances, [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) has more false positives at longer distances. Real Hi-C data normalized using KR and MA normalizations. [R package](https://bitbucket.org/nadhir/find/downloads/). <details>
    <summary>Paper</summary>
    Mohamed Nadhir, Djekidel, Yang Chen, and Michael Q. Zhang. “FIND: DifFerential Chromatin INteractions Detection Using a Spatial Poisson Process.” Genome Research, February 12, 2018. https://doi.org/10.1101/gr.212241.116.
</details>

- <a name="diffloop">[diffloop](https://bioconductor.org/packages/diffloop/)</a> - Differential analysis of chromatin loops (ChIA-PET). [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html) framework. <details>
    <summary>Paper</summary>
    Lareau, Caleb A., and Martin J. Aryee. "Diffloop: A Computational Framework for Identifying and Analyzing Differential DNA Loops from Sequencing Data"  https://doi.org/10.1093/bioinformatics/btx623 Bioinformatics (Oxford, England), September 29, 2017. 
</details>

- <a name="ap">[AP](https://github.com/XiaoTaoWang/TADLib)</a> - aggregation preference - parameter, to quantify TAD heterogeneity. Call significant interactions within a TAD, cluster with DBSCAN, calculate weighted interaction density within each cluster, average. AP measures are reproducible. Comparison of TADs in Gm12878 and IMR90 - stable TADs change their aggregation preference, these changes correlate with LINEs, Lamin B1 signal. Can detect structural changes (block split) in TADs. <details>
    <summary>Paper</summary>
    Wang, X.-T., Dong, P.-F., Zhang, H.-Y., and Peng, C. (2015). "Structural heterogeneity and functional diversity of topologically associating domains in mammalian genomes.](https://academic.oup.com/nar/article/43/15/7237/2414371)" Nucleic Acids Research
</details>

- <a name="diffhic">[diffHiC](https://bioconductor.org/packages/diffHic/)</a> - Differential contacts using the full pipeline for Hi-C data. Explanation of the technology, binning. MA normalization, [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)-based. Comparison with [HOMER](#homer). [Documentation](https://bioconductor.org/packages/release/bioc/vignettes/diffHic/inst/doc/diffHic.pdf). <details>
    <summary>Paper</summary>
    Lun, Aaron T. L., and Gordon K. Smyth. "DiffHic: A Bioconductor Package to Detect Differential Genomic Interactions in Hi-C Data"  https://doi.org/10.1186/s12859-015-0683-0 BMC Bioinformatics 16 (2015)
</details>

- <a name="hiccupsdiff">[HiCCUPS Diff](https://github.com/aidenlab/juicer/wiki/HiCCUPSDiff)</a> - differential loop analysis. Input - two .hic files and loop lists; output - lists of differential loops. 

## TAD callers

- <a name="supertld">[SuperTLD](https://github.com/deepomicslab/SuperTLD)</a> - hierarchical TAD-like domain (TLD) detection from RNA-DNA associated interactions (RAIs, technologies like iMARGI, GRID-seq, ChAR-seq, RADICL-seq). Adapts [SuperTAD](#supertad) by incorporating a Bayesian correction during dynamic programming and iterative inferring hierarchy. Incorporates imputation of missing interaction frequencies using a negative binomial model (inspired by [SAVER](https://github.com/mohuangx/SAVER)). Applied on four RAI datasets, compares TLDs with TADs (Overlapping Ratio (OR) measuring similarity between two hierarchies, normalized mutual information (NMI) measuring similarity of two disjoint partitions (metrics from [SuperTAD](#supertad) paper, four types of TAD relationships (matched, merged, split, shifted, adapted from [SpectralTAD](#spectraltad)), also Pearson correlation and distance decay). TLDs and TADs are moderately similar. RAI-constructed similarity matrices and Hi-C maps are highly correlated, TADs and TLDs are moderately similar, more novel boundaries in TLDs. Common boundaries are associated with cytoplasm and cytosol terms, TLD boundaries - with different terms. <details>
  <summary>Paper</summary>
  Zhang, Yu Wei, Lingxi Chen, and Shuai Cheng Li. "Detecting TAD-like domains from RNA-associated interactions." Nucleic Acids Research (25 May 2022). https://doi.org/10.1093/nar/gkac422
</details>

- <a name="supertad">[SuperTAD](https://zenodo.org/record/4314123)</a> - hierarchical TAD callers, efficient algorithms for computing the coding tree of a contact map, structure information theory, dynamic programming, polynomial-time solvable. Compared with seven hierarchical TAD callers (but not [SpectralTAD](#spectraltad)!), better epigenomic enrichment, agreement in calls from raw and KR matrices, better cross-resolution agreement by size, width. Command line, C++ implementation. <details>
  <summary>Paper</summary>
  Zhang, Yu Wei, Meng Bo Wang, and Shuai Cheng Li. “SuperTAD: Robust Detection of Hierarchical Topologically Associated Domains with Optimized Structural Information.” Genome Biology 22, no. 1 (December 2021): 45. https://doi.org/10.1186/s13059-020-02234-6.
</details>

- <a name="grinch">[GRiNCH](https://roy-lab.github.io/grinch/)</a> - TAD detection algorithm using Graph Regularized Nonnegative matrix factorization. Graph captures distance dependence. Also smoothes/imputes Hi-C matrices. Compared with Directionality Score, Insulation Index, [rGMAP](https://github.com/wbaopaul/rGMAP), [Armatus](#armatus), [HiCSeg](#hicseg), [TopDom](#topdom) using several metrics - Davies-Bouldin index, Delta Contact Counts, Rand Index, Mutual Information. GRiNCH detects consistent similarly-sized TADs, robust to different resolutions, boundaries show robust enrichment in known markers of TAD boundaries (TFBSs and histone), detects consistent [Fit-Hi-C](#fit-hi-c) significant interactions (area under precision-recall curve). Applied to mouse neural development and pluripotency reprogramming to confirm known and discover new boundary regulators. Applied to SPRITE and HiChIP data. [Visualization tutorial](https://github.com/Roy-lab/grinch/tree/master/visualization). <details>
    <summary>Paper</summary>
    Lee, Da-Inn, and Sushmita Roy. "Graph-Regularized Matrix Factorization for Reliable Detection of Topological Units from High-Throughput Chromosome Conformation Capture Datasets"  https://doi.org/10.1186/s13059-021-02378-z Genome Biology, 25 May 2021
</details>

- <a name="hickey">[HiCKey](https://github.com/YingruWuGit/HiCKey)</a> - hierarchical TAD caller, comparison of TADs across samples. A generalized likelihood-ratio (GRL) test for detecting change-points in an interaction matrix that follows negative binomial distribution (Methods). Bottom-up approach to detect hierarchy. Tested on Forcato simulated data with nested TADs, TPR/FPR/difference/Fowlkes-Mallows index to estimate performance. Applied to seven cell lines. TAD hierarchy is up to four levels. Compared with [TADtree](#tadtree), [3D-NetMod](#3d-netmod), [IC-Finder](#ic-finder), [HiCSeg](#hicseg). Colocalization within a 2-bin distance. Input - normalized, distance effect removed matrix in sparse text format, output - TAD start coordinate, hierarchy level, p-value of the changepoint. C++ implementation. Did not compare with [SpectralTAD](https://bioconductor.org/packages/SpectralTAD/) hierarchical caller.<details>
    <summary>Paper</summary>
    Xing, Haipeng. "Deciphering Hierarchical Organization of Topologically Associated Domains through Change-Point Testing"  https://doi.org/10.1186/s12859-021-04113-8 BMC Bioinformatics, April 10, 2021
</details>

- <a name="tadbd">[TADBD](https://github.com/bioinfo-lab/TADBD/)</a> - TAD caller using a multi-scale Haar diagonal template (sum of on-diagonal squares minus the sum of off-diagonal squares). Compared with [HiCDB](#hicdb), [IC-Finder](#ic-finder), [EAST](#east) (also using Haar features), [TopDom](#topdom), [HiCseg](#hicseg) using simulated (Forcato) and experimental (K562 and IMR90) data. ICE-normalized data. MCC, Jaccard. [R package](https://rdrr.io/github/bioinfo-lab/TADBD/src/R/TADBD.R). <details>
    <summary>Paper</summary>
    Lyu, Hongqiang, Lin Li, Zhifang Wu, Tian Wang, Jiguang Zheng, and Hongda Wang. "TADBD: A Sensitive and Fast Method for Detection of Typologically Associated Domain Boundaries"  https://doi.org/10.2144/btn-2019-0165 BioTechniques, April 7, 2020
</details>

- <a name="bhi-cect">[BHi-Cect](https://github.com/princeps091-binf/BHi-Cect)</a> - identification of the full hierarchy of chromosomal interactions (TADs). Spectral clustering starting from the whole chromosome, detecting nested BHi-Cect Partition Trees (BPTs), partitioned in non-contiguous and interwoven enclaves, inspired by fractal globule idea. Variation of information to test the agreement between two clustering results, overlap-based metrics to test correspondence with TADs. Correspondence analysis of enclaves association with TF content. Gene enrichment. Different enclaves show different epigenomic and gene expression signatures, bottom enclaves are most crisply defined. Resolution affects what enclave size can be detected. <details>
    <summary>Paper</summary>
    Kumar, Vipin, Simon Leclerc, and Yuichi Taniguchi. "BHi-Cect: A Top-down Algorithm for Identifying the Multi-Scale Hierarchical Structure of Chromosomes"  https://doi.org/10.1093/nar/gkaa004 Nucleic Acids Research 48, no. 5 (March 18, 2020)
</details>

- <a name="tadpole">[TADpole](https://github.com/3DGenomes/TADpole)</a> - hierarchical TAD boundary caller. Preprocessing by filtering sparse rows, transforming the matrix into its Pearson correlation coefficient matrix, running PCA on it and retaining 200 PCs, transforming into a Euclidean distance matrix, clustering using the Constrained Incremental Sums of Squares clustering (`rioja::chclust(, coniss)`), estimating significance, Calinski-Harabasz index to estimate the optimal number of clusters (chromatin subdivisions). Benchmarking using Zufferey 2018 datasets, mouse limb bud development with genomic inversions from Kraft 2019. Resolution, normalization, sequencing depth. Metrics: the Overlap Score, the Measure of Concordance, all from Zufferey 2018. Enrichment in epigenomic marks. DiffT metric for differential analysis (on binarized TAD/non-TAD matrices). Compared with 22 TAD callers, including hierarchical ([CaTCH](#catch), [rGMAP](https://github.com/wbaopaul/rGMAP), [Matryoshka](https://github.com/COMBINE-lab/matryoshka), [PSYCHIC](https://github.com/dhkron/PSYCHIC)). <details>
    <summary>Paper</summary>
    Soler-Vila, Paula, Pol Cuscó Pons, Irene Farabella, Marco Di Stefano, and Marc A. Marti-Renom. "Hierarchical Chromatin Organization Detected by TADpole"  https://doi.org/10.1101/698720 Preprint. Bioinformatics, July 11, 2019. 
</details>

- <a name="hicdb">[HiCDB](https://github.com/ChenFengling/RHiCDB)</a> - TAD boundary detection using local relative insulation (LRI) metric, improved stability, less parameter tuning, cross-resolution, differential boundary detection, lower computations, visualization. Review of previous methods, directionality index, insulation score. Math of LRI. GSEA-like enrichment in genome annotations (CTCF). Differential boundary detection using the intersection of extended boundaries. Compared with [Armatus](#armatus), [DomainCaller](#domaincaller), [HiCseg](#hicseg), [IC-Finder](#ic-finder), Insulation, [TopDom](#topdom) on 40kb datasets. Accurately detects smaller-scale boundaries. Differential TADs are enriched in cell-type-specific genes. <details>
    <summary>Paper</summary>
    Chen, Fengling, Guipeng Li, Michael Q. Zhang, and Yang Chen. "HiCDB: A Sensitive and Robust Method for Detecting Contact Domain Boundaries"  https://doi.org/10.1093/nar/gky789 Nucleic Acids Research 46, no. 21 (November 30, 2018)
</details>

- <a name="ontad">[OnTAD](https://github.com/anlin00007/OnTAD)</a> - hierarchical TAD caller, Optimal Nested TAD caller. Sliding window, adaptive local minimum search algorithm, similar to [TOPDOM](#topdom). Other hierarchical callers - [TADtree](#tadtree), [rGMAP](https://github.com/wbaopaul/rGMAP), [Arrowhead](#arrowhead), [3D-Net](#3d-net), [IC-Finder](#ic-finder). Performance comparison with [DomainCaller](#domaincaller), [rGMAP](https://github.com/wbaopaul/rGMAP), [Arrowhead](#arrowhead), [TADtree](#tadtree). Stronger enrichment of CTCF and two cohesin proteins RAD21 and SMC3. C++ implementation. [OnTAD for coolers](https://github.com/cchlanger/cooler_ontad) - a Python wrapper to work with `.cool` files. <details>
    <summary>Paper</summary>
    An, Lin, Tao Yang, Jiahao Yang, Johannes Nuebler, Qunhua Li, and Yu Zhang. "Hierarchical Domain Structure Reveals the Divergence of Activity among TADs and Boundaries"  https://doi.org/10.1101/361147 July 3, 2018. - Intro about TADs, Dixon's directionality index, Insulation score. Other hierarchical callers - TADtree, rGMAP, Arrowhead, 3D-Net, IC-Finder. Limitations of current callers - ad hoc thresholds, sensitivity to sequencing depth and mapping resolution, long running time and large memory usage, insufficient performance evaluation. Boundaries are asymmetric - some have more contacts with other boundaries, support for asymmetric loop extrusion model. Performance comparison with DomainCaller, rGMAP, Arrowhead, TADtree. Stronger enrichment of CTCF and two cohesin proteins RAD21 and SMC3. TAD-adjR^2 metric quantifying the proportion of variance in the contact frequencies explained by TAD boundaries. Reproducibility of TAD boundaries - Jaccard index, tested at different sequencing depths and resolutions. Boundaries of hierarchical TADs are more active - more CTCF, epigenomic features, TFBSs expressed genes. Super-boundaries - shared by 5 or more TADs, highly active. Rao-Huntley 2014 GM12878 data. Distance correction - subtracting the mean counts at each distance.
</details>

- <a name="3d-netmod">[3D-NetMod](https://bitbucket.org/creminslab/3dnetmod_method_v1.0_10_06_17)</a> - hierarchical, nested, partially overlapping TAD detection using graph theory. Community detection method based on the maximization of network modularity, Louvain-like locally greedy algorithm, repeated several (20) times to avoid local maxima, then getting consensus. Tuning parameters are estimated over a sequence search. Benchmarked against [TADtree](#tadtree), directionality index, [Arrowhead](#arrowhead). ICE-normalized data brain data from Geschwind (human data) and Jiang (mouse data) studies. Computationally intensive. Python implementation. <details>
    <summary>Paper</summary>
    Norton, Heidi K., Daniel J. Emerson, Harvey Huang, Jesi Kim, Katelyn R. Titus, Shi Gu, Danielle S. Bassett, and Jennifer E. Phillips-Cremins. “Detecting Hierarchical Genome Folding with Network Modularity.” Nature Methods 15, no. 2 (February 2018): 119–22. https://doi.org/10.1038/nmeth.4560.
</details>

- <a name="dedoc">[deDoc](https://github.com/yinxc/structural-information-minimisation)</a> - TAD detection minimizing structural entropy of the Hi-C graph (structural information theory). Detects optimal resolution (= minimal entropy). Pooled 10 single-cell Hi-C analysis. Intro about TADs, a brief description of TAD callers, including hierarchical. Works best on raw, non-normalized data, highly robust to sparsity (0.1% of the original data sufficient). Compared with five TAD callers ([Armatus](#armatus), [TADtree](#tadtree), [Arrowhead](#arrowhead), [MrTADFinder](https://github.com/gersteinlab/MrTADFinder), [DomainCaller (DI)](#domaincaller)), and a classical graph modularity detection algorithm. Enrichment in CTCF, housekeeping genes, H3K4me3, H4K20me1, H3K36me3. Other benchmarks - weighted similarity, number, length of TADs. Detects hierarchy over different passes. Java implementation (won't run on Mac). <details>
    <summary>Paper</summary>
    Li, Angsheng, Xianchen Yin, Bingxiang Xu, Danyang Wang, Jimin Han, Yi Wei, Yun Deng, Ying Xiong, and Zhihua Zhang. “Decoding Topologically Associating Domains with Ultra-Low Resolution Hi-C Data by Graph Structural Entropy.” Nature Communications 9, no. 1 (15 2018): 3265. https://doi.org/10.1038/s41467-018-05691-7.
</details>

- <a name="catch">[CaTCH](https://github.com/zhanyinx/CaTCH_R)</a> - identification of hierarchical TAD structure. Reciprocal insulation (RI) index. Benchmarked against Dixon's TADs (diTADs). CTCF enrichment as a benchmark, enrichment of TADs in differentially expressed genes. <details>
    <summary>Paper</summary>
    Zhan, Yinxiu, Luca Mariani, Iros Barozzi, Edda G. Schulz, Nils Blüthgen, Michael Stadler, Guido Tiana, and Luca Giorgetti. "Reciprocal Insulation Analysis of Hi-C Data Shows That TADs Represent a Functionally but Not Structurally Privileged Scale in the Hierarchical Folding of Chromosomes"  https://doi.org/10.1101/gr.212803.116 Genome Research 27, no. 3 (2017)
</details>

- <a name="hitad">[HiTAD](https://github.com/XiaoTaoWang/TADLib)</a> - hierarchical TAD identification, different resolutions, correlation with chromosomal compartments, replication timing, gene expression. Adaptive directionality index approach. Data sources, methods for comparing TAD boundaries, reproducibility. H3K4me3 enriched and H3K4me1 depleted at boundaries. TAD boundaries (but not sub-TADs) separate replication timing, A/B compartments, gene expression. <details>
    <summary>Paper</summary>
    Wang, Xiao-Tao, Wang Cui, and Cheng Peng. "HiTAD: Detecting the Structural and Functional Hierarchies of Topologically Associating Domains from Chromatin Interactions"  https://doi.org/10.1093/nar/gkx735 Nucleic Acids Research 45, no. 19 (November 2, 2017)
</details>

- <a name="clustertad">[ClusterTAD](https://github.com/BDM-Lab/ClusterTAD)</a> - A clustering method for identifying topologically associated domains (TADs) from Hi-C data. <details>
    <summary>Paper</summary>
    Oluwadare, Oluwatosin, and Jianlin Cheng. "ClusterTAD: An Unsupervised Machine Learning Approach to Detecting Topologically Associated Domains of Chromosomes from Hi-C Data" https://doi.org/10.1186/s12859-017-1931-2 BMC Bioinformatics 18, no. 1 (November 14, 2017)
</details>

- <a name="ic-finder">[IC-Finder](https://academic.oup.com/nar/article/45/10/e81/2958551)</a> - Segmentations of HiC maps into hierarchical interaction compartments. Code not available. <details>
    <summary>Paper</summary>
    Noelle Haddad, Cedric Vaillant, Daniel Jost. "IC-Finder: inferring robustly the hierarchical organization of chromatin folding" https://doi.org/10.1093/nar/gkx036 
</details>

- <a name="east">[EAST](https://github.com/ucrbioinfo/EAST)</a> - Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps, Haar-like features (rectangles on images) and a function that quantifies TAD properties: frequency within is high, outside - low, boundaries must be strong. Objective - finding a set of contiguous non-overlapping domains maximizing the function. Restricted by the maximum length of TADs. Boundaries are enriched in CTCF, RNP PolII, H3K4me3, H3K27ac. <details>
    <summary>Paper</summary>
    Abbas Roayaei Ardakany, Stefano Lonardi, and Marc Herbstritt, "Efficient and Accurate Detection of Topologically Associating Domains from Contact Maps"  https://doi.org/10.4230/LIPIcs.WABI.2017.22 (Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik GmbH, Wadern/Saarbruecken, Germany, 2017)
</details>

- <a name="lavaburst">[lavaburst](https://github.com/nvictus/lavaburst)</a> - extension of the [Armatus](#armatus) algorithm introduced in [Filippova et al., “Identification of Alternative Topological Domains in Chromatin”](https://doi.org/10.1186/1748-7188-9-14). [Documentation](https://nezar-compbio.github.io/lavaburst/index.html). <details>
    <summary>Paper</summary>
    Schwarzer, Wibke, Nezar Abdennur, Anton Goloborodko, Aleksandra Pekowska, Geoffrey Fudenberg, Yann Loe-Mie, Nuno A. Fonseca, et al. “Two Independent Modes of Chromatin Organization Revealed by Cohesin Removal.” Nature 551, no. 7678 (02 2017): 51–56. https://doi.org/10.1038/nature24281.
</details>

- <a name="arboretum-hi-c">[Arboretum-Hi-C](https://bitbucket.org/roygroup/arboretum-hic)</a> - a multitask spectral clustering method to identify differences in genomic architecture. Intro about the 3D genome organization, TAD differences, and conservation. Assessment of different clustering approaches using different distance measures, as well as raw contacts. Judging clustering quality by enrichment in genomic regulatory signals (Histone marks, LADs, early vs. late replication timing, TFs like POLII, TAF, TBP, CTCF, P300, CMYC, cohesin components, LADs, replication timing, SINE, LINE, LTR) and by numerical methods (Davies-Bouldin index, silhouette score, others). Although spectral clustering on contact counts performed best, spectral + Spearman correlation was chosen. Comparing cell types identifies biologically relevant differences as quantified by enrichment. Peak counts or average signal within regions were used for enrichment. [Data](https://zenodo.org/record/49767). <details>
    <summary>Paper</summary>
    Fotuhi Siahpirani, Alireza, Ferhat Ay, and Sushmita Roy. "A Multi-Task Graph-Clustering Approach for Chromosome Conformation Capture Data Sets Identifies Conserved Modules of Chromosomal Interactions"  https://doi.org/10.1186/s13059-016-0962-8 Genome Biology 17, no. 1 (December 2016). 
</details>

- <a name="tadtree">[TADtree](https://github.com/raphael-group/TADtree)</a> - Hierarchical (nested) TAD identification. Two ways of TAD definition: 1D and 2D. Normalization by distance. Enrichment over the background. [Documentation](http://compbio.cs.brown.edu/projects/tadtree/) <details>
    <summary>Paper</summary>
    Weinreb, Caleb, and Benjamin J. Raphael. "Identification of Hierarchical Chromatin Domains"  https://doi.org/10.1093/bioinformatics/btv485 Bioinformatics (Oxford, England) 32, no. 11 (June 1, 2016)
</details>

- <a name="topdom">[TopDom](https://github.com/HenrikBengtsson/TopDom)</a> - An efficient and Deterministic Method for identifying Topological Domains in Genomes, Method is based on the general observation that within-TAD interactions are stronger than between-TAD. binSignal value as the average of nearby contact frequency, fitting a curve, finding local minima, test them for significance. Fast, takes linear time. Detects similar domains to HiCseq and Dixon's directionality index. Found expected enrichment in CTCF, histone marks. Housekeeping genes and overall gene density are close to TAD boundaries, differentially expressed genes are not. [R package](https://cran.r-project.org/web/packages/TopDom/index.html). <details>
    <summary>Paper</summary>
    Shin, Hanjun, Yi Shi, Chao Dai, Harianto Tjong, Ke Gong, Frank Alber, and Xianghong Jasmine Zhou. "TopDom: An Efficient and Deterministic Method for Identifying Topological Domains in Genomes"  https://doi.org/10.1093/nar/gkv1505 Nucleic Acids Research 44, no. 7 (April 20, 2016)
</details>

- <a name="tadtool">[TADtool](https://github.com/vaquerizaslab/tadtool)</a> - command-line wrapper for directionality index and insulation score TAD callers. <details>
    <summary>Paper</summary>
    Kruse, Kai, Clemens B. Hug, Benjamín Hernández-Rodríguez, and Juan M. Vaquerizas. "TADtool: Visual Parameter Identification for TAD-Calling Algorithms"  https://doi.org/10.1093/bioinformatics/btw368 Bioinformatics (Oxford, England) 32, no. 20 (15 2016)
</details>

- <a name="armatus">[Armatus](https://github.com/kingsfordgroup/armatus)</a> - TAD detection at different resolutions. Dynamic programming method. <details>
    <summary>Paper</summary>
    Filippova, Darya, Rob Patro, Geet Duggal, and Carl Kingsford. "Identification of Alternative Topological Domains in Chromatin" Algorithms for Molecular Biology 9, no. 1 (2014) https://doi.org/10.1186/1748-7188-9-14
</details>

- <a name="hicseg">[HiCseg](https://cran.r-project.org/web/packages/HiCseg/index.html)</a> - TAD detection by maximization of likelihood based block-wise segmentation model. 2D segmentation rephrased as 1D segmentation - not contours, but borders. Statistical framework, solved with dynamic programming. Dixon data as gold standard. [Hausdorff distance](https://en.wikipedia.org/wiki/Hausdorff_distance) to compare segmentation quality. Parameters (from [TopDom](#topdom) paper): nb_change_max = 500, distrib = 'G' and model = 'Dplus'. <details>
    <summary>Paper</summary>
    Lévy-Leduc, Celine, M. Delattre, T. Mary-Huard, and S. Robin. "Two-Dimensional Segmentation for Analyzing Hi-C Data"  https://doi.org/10.1093/bioinformatics/btu443 Bioinformatics (Oxford, England) 30, no. 17 (September 1, 2014)
</details>

- <a name="domaincaller">[domaincaller](https://github.com/XiaoTaoWang/domaincaller)</a> - A Python implementation of the original DI domain caller.

- <a name="arrowhead">[Arrowhead](https://github.com/aidenlab/juicer/wiki/Arrowhead)</a> - contact domain (TAD) detection using Arrowhead transformation. Described in Section IV.a of the Extended Experimental Procedures of [Rao, Huntley et al. Cell 2014](https://www.cell.com/cms/10.1016/j.cell.2014.11.021/attachment/d3c6dcd8-c799-4f68-bbe4-201be54960b5/mmc1.pdf)

### TAD detection, benchmarking

- [TADMaster](http://biomlearn.uccs.edu/TADMaster/) - a tool/webserver to evaluate the concordance of results of TAD callers. Multiple comparison metrics (number, size, Measure of Concordance (MoC) considering tolerance boundary (flank)) Two run modes - one accepts BED files and compares them, another (TADMaster Plus) processes Hi-C matrices (hic, h5, cool, sparse, full matrices), normalizes them (VC, ICE, KR, more), detects TADs (DI, Insulation Score, HiCseq, SpectralTAD, 12 total), and compares them. Visualization, clustering, [example output](http://biomlearn.uccs.edu/TADMaster/visualize_example/415/). Docker version, [GitHub](https://github.com/OluwadareLab/TADMaster). <details>
  <summary>Paper</summary>
  Higgins, Sean, Victor Akpokiro, Allen Westcott, and Oluwatosin Oluwadare. “TADMaster: A Comprehensive Web-Based Tool for the Analysis of Topologically Associated Domains.” BMC Bioinformatics 23, no. 1 (November 4, 2022): 463. https://doi.org/10.1186/s12859-022-05020-2.
</details>

- Sefer, Emre. “[A Comparison of Topologically Associating Domain Callers over Mammals at High Resolution](https://doi.org/10.1186/s12859-022-04674-2).” BMC Bioinformatics 23, no. 1 (December 2022) - Benchmarking of 27 TAD callers (Table 1). Brief description of each caller. Classified into 3 categories: feature-based, clustering, graph-partitioning methods. Systematic evaluation of TAD number, size, enrichment/depletion in functional annotations, 8 metrics. Distinguish corner-dot loops and TADs withour corner dots and evaluate performance on detecting each. Using mESCs Micro-C (HiCNorm normalized), simulated data, testing resolution, sequencing depth. Feature-based callers generally perform better, but performance depends on the type of test. SpectralTAD performs well. Code for all tests at [GitHub](https://github.com/seferlab/TADcomparison).

- [Brief description of 22 TAD calling methods](https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-018-1596-9/MediaObjects/13059_2018_1596_MOESM1_ESM.pdf). Source: [Zufferey et al., “Comparison of Computational Methods for the Identification of Topologically Associating Domains.”](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1596-9#Bib1)

- Zufferey, Marie, Daniele Tavernari, Elisa Oricchio, and Giovanni Ciriello. "[Comparison of Computational Methods for the Identification of Topologically Associating Domains](https://doi.org/10.1186/s13059-018-1596-9)" Genome Biology 19, no. 1 (10 2018) - Comparison of 22 TAD callers across different conditions. Callers are classified as linear score-based, statistical model-based, clustering, graph theory. Table 1 and Additional file 1 summarizes each caller. The effect of data resolution, normalization, hierarchy. Test on Rao 2014 data, chromosome 6. ICE or LGF (local genomic feature) normalization. The measure of Concordance (MoC) to compare TADs. CTCF/cohesin as a measure of biological significance. TopDom, HiCseg, CaTCH, CHDF are the top performers. [R scripts, including for calculation MoC](https://github.com/CSOgroup/TAD-benchmarking-scripts)

- Dali, Rola, and Mathieu Blanchette. "[A Critical Assessment of Topologically Associating Domain Prediction Tools](https://doi.org/10.1093/nar/gkx145)" Nucleic Acids Research 45, no. 6 (April 7, 2017) - TAD definition, tools. Meta-TADs, hierarchy, overlapping TADs.[HiCPlotter](#hicplotter) for visualization. Manual annotation as a gold standard. Sequencing depth and resolution affect things. Code, manual annotations

- Forcato, Mattia, Chiara Nicoletti, Koustav Pal, Carmen Maria Livi, Francesco Ferrari, and Silvio Bicciato. "[Comparison of Computational Methods for Hi-C Data Analysis](https://doi.org/10.1038/nmeth.4325)" Nature Methods, June 12, 2017 - Hi-C processing and TAD calling tools benchmarking, Table 1, simulated (Lun and Smyth method) and real data. Notes about pluses and minuses of each tool. TAD reproducibility is higher than chromatin interactions, increases with a larger number of reads. Consistent enrichment of TAD boundaries in CTCF, irrespectively of TAD caller. Hi-C replication is poor, just a bit more than random. Supplementary table 2 - technical details about each program, Supplementary Note 1 - Hi-C preprocessing tools, Supplementary Note 2 - TAD callers. [Supplementary note 3](https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf) - how to simulate Hi-C data. [Supplementary note 6](https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf) - how to install tools. [Tools for TAD comparison, and simulated matrices](https://bitbucket.org/mforcato/hictoolscompare.git).

- Olivares-Chauvet, Pedro, Zohar Mukamel, Aviezer Lifshitz, Omer Schwartzman, Noa Oded Elkayam, Yaniv Lubling, Gintaras Deikus, Robert P. Sebra, and Amos Tanay. "[Capturing Pairwise and Multi-Way Chromosomal Conformations Using Chromosomal Walks](https://doi.org/10.1038/nature20158)" Nature 540, no. 7632 (November 30, 2016) - TADs organize chromosomal territories. Active and inactive TAD properties. Methods: Good mathematical description of insulation score calculations. Filter TADs smaller than 250kb. Inter-chromosomal contacts are rare, \~7-10%. Concatemers (more than two contacts) are unlikely.

- Rocha, Pedro P., Ramya Raviram, Richard Bonneau, and Jane A. Skok. "[Breaking TADs: Insights into Hierarchical Genome Organization](https://doi.org/10.2217/epi.15.25)" Epigenomics 7, no. 4 (2015) - Textbook overview of TADs in 3 pages with key references. 3D organization discovery using FISH, 3C, Hi-C. Discovery of A/B compartments (euchromatin, heterochromatin), TADs as regulatory units conserved even across syntenic regions in different organisms. TADs coordinate gene expression. TAD boundaries are not created equal. Examples of changes of TAD boundaries (Hox gene cluster, ES differentiation). Hierarchy of TADs.

- Crane, Emily, Qian Bian, Rachel Patton McCord, Bryan R. Lajoie, Bayly S. Wheeler, Edward J. Ralston, Satoru Uzawa, Job Dekker, and Barbara J. Meyer. "[Condensin-Driven Remodelling of X Chromosome Topology during Dosage Compensation](https://doi.org/10.1038/nature14450)" Nature 523, no. 7559 (July 9, 2015). - Insulation Score to define TADs - sliding square along the diagonal, aggregating signal within it. This aggregated score is normalized and binned into TADs, boundaries. See [Methods and implementation](https://github.com/dekkerlab/crane-nature-2015). [matrix2insulation.pl](https://github.com/dekkerlab/cworld-dekker/tree/master/scripts/perl/matrix2insulation.pl), Parameters: -is 480000 -ids 320000 -im iqrMean -nt 0 -ss 160000 -yb 1.5 -nt 0 -bmoe 0. 

- "Hierarchical Regulatory Domain Inference from Hi-C Data" - presentation by Bartek Wilczyński about TAD detection, existing algorithms, new SHERPA and OPPA methods. [Video](https://simons.berkeley.edu/talks/bartek-wilczynski-03-10-16), [PDF](https://simons.berkeley.edu/sites/default/files/docs/4588/2016-03-10-simons-institute-wilczynski.pdf), [Web site](http://regulomics.mimuw.edu.pl/wp/), [GitHub](https://github.com/regulomics/) - SHERPA and OPPA code there.

### Architectural stripes

- <a name="stripenn">[Stripenn](https://github.com/vahedilab/stripenn)</a> - a computer vision-based method for architectural stripes detection using Canny edge detection.Scores stripes by median p-value and stripiness based on the continuity of interaction signal. Input - .cool files, optionally normalized. Output - coordinates and scores of the predicted stripes. Applicable to Hi-C, HiChIP, Micro-C data. Introduction to the biology of architectural stripes, review of previous methods (Zebra from Vian et al. 2018, [domainClassifyR](https://github.com/ChristopherBarrington/domainClassifyR), [CHESS](#chess) for comparing 3D domains and stripes). Analysis of stripes from B and T lymphocytes identifies stripe anchors enriched in the transcriptionally active compartments, architectural proteins mediating loop extrusion. Strips are strongly conserved, correspond to TAD boundaries, subtle changes are associated with transcriptional output. Python, three functions (`compute`, `score`, `seeimage`). [Video 16m](https://www.youtube.com/watch?v=_TEomDYv2vk). <details>
    <summary>Paper</summary>
    Yoon, Sora, and Golnaz Vahedi. "Stripenn Detects Architectural Stripes from Chromatin Conformation Data Using Computer Vision"  https://doi.org/10.1101/2021.04.16.440239 Preprint. Bioinformatics, April 18, 2021
</details>

- [domainClassifyR](https://github.com/ChristopherBarrington/domainClassifyR) - detection and classification of 'Stripes' structures at TAD boundaries. Associated with both poised and active chromatin landscapes and transcription is not a key determinant of their structure. Stripe formation is linked to the functional state of the cell through cohesin loading at lineage-specific enhancers and developmental control of CTCF binding site occupancy. Comparison of pluripotent mouse embryonic stem cells and lineage-committed neural cells and characterizing emergent, lost, and mainained stripes. Raw data at [SRA668328](https://www.ncbi.nlm.nih.gov/sra/?term=SRA668328). <details>
    <summary>Paper</summary>
    Barrington, Christopher, Dimitra Georgopoulou, Dubravka Pezic, Wazeer Varsally, Javier Herrero, and Suzana Hadjur. "Enhancer accessibility and CTCF occupancy underlie asymmetric TAD architecture and cell type specific genome topology." Nature communications 10, no. 1 (2019): 1-14. https://doi.org/10.1038/s41467-019-10725-9
</details>

- <a name="stripecaller">[StripeCaller](https://github.com/XiaoTaoWang/StripeCaller)</a> - A toolkit for analyzing architectural stripes. Architectural stripes, created by extensive loading of cohesin near CTCF anchors, with Nipbl and Rad21 help. Little overlap between B cells and ESCs. Architectural stripes are sites for tumor-inducing TOP2beta DNA breaks. ATP is required for loop extrusion, cohesin translocation, but not required for maintenance, Replication of transcription is not important for loop extrusion. Zebra algorithm for detecting architectural stripes, image analysis, math in Methods. Human lymphoblastoid cells, mouse ESCs, mouse B-cells activated with LPS, CH12 B lymphoma cells, wild-type, treated with hydroxyurea (blocks DNA replication), flavopiridol (blocks transcription, PolII elongation), oligomycin (blocks ATP). Hi-C, ChIA-pet, ChIP-seq, ATAC-seq, and more [Data1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE82144), [Data2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98119). <details>
    <summary>Paper</summary>
    Vian, Laura, Aleksandra Pękowska, Suhas S.P. Rao, Kyong-Rim Kieffer-Kwon, Seolkyoung Jung, Laura Baranello, Su-Chen Huang, et al. "The Energetics and Physiological Impact of Cohesin Extrusion"  https://doi.org/10.1016/j.cell.2018.03.072 Cell 173, no. 5 (May 2018)
</details>


### Differential, timecourse TAD analysis

- [DiffDomain](https://github.com/Tian-Dechao/diffDomain) - differential TAD analysis. Uses random matrix theory (the Tracy–Widom distribution) to identify structurally reorganized TADs (six subtypes: strength-change, loss, split, merge, zoom, and complex). Outperforms other methods (TADCompare, HiCcompare, DiffGR, DiffTAD, HiC-DC+, TADsplimer) in terms of FPR, TPR, proportions and subtypes of reorganized TADs. Applicable to scHi-C (pseudobulk, 100 cells minimum, imputed with scHiCluster). Input - two Hi-C matrices and TADs called in the first (e.g., Arrowhead on KR-normalized matrix). Computes a log2 difference matrix from two KR-normalized matrices, then distance/depth-normalizes off-diagonals, normalizes by SQRT(N) - white noise, a generalized Wigner D matrix, the largest eigenvalue of which converges to 2, hypothesis testing is the eigenvalue is larger than 2. <details>
  <summary>Paper</summary>
  Hua, Dunming, Ming Gu, Yanyi Du, Li Qi, Xiangjun Du, Zhidong Bai, Xiaopeng Zhu, and Dechao Tian. "DiffDomain enables identification of structurally reorganized topologically associating domains." In International Conference on Research in Computational Molecular Biology, pp. 302-303. Springer, Cham, 2022. https://doi.org/10.1007/978-3-031-04749-7_20
</details>

- <a name="intertads">[InterTADs](https://github.com/nikopech/InterTADs)</a> - TAD-centric integration of multi-omics (SNPs, gene expression, methylation) data by genomic coordinateы, proper data scaling, and differential analysis. Differential TAD analysis by interaction strength changes, the number of other differential events. [Documentation](https://github.com/nikopech/InterTADs/wiki/Usage). <details>
    <summary>Paper</summary>
    - Tsagiopoulou, Maria, Nikolaos Pechlivanis, and Fotis Psomopoulos. "InterTADs: Integration of Multi-Omics Data on Topological Associated Domains"  https://doi.org/10.21203/rs.3.rs-54194/v1 Preprint. In Review, August 12, 2020
</details>

- <a name="tadcompare">[TADcompare](https://github.com/dozmorovlab/TADCompare)</a> - R package for differential and time-course TAD boundary analysis. Uses [SpectralTAD](https://bioconductor.org/packages/SpectralTAD/) score - spectral decomposition of Hi-C matrices - to statistically detect five types of differential TAD boundaries: merge, split, complex, shifted, strength change. In the time-course analysis, detects six types of boundary score changes: highly common, early appearing, late appearing, early disappearing, late disappearing, and dynamic TAD boundaries. Returns genomic coordinated and types of TAD boundary changes in BED format. [Documentation](https://dozmorovlab.github.io/TADCompare/), [Bioconductor Package](https://bioconductor.org/packages/TADCompare/) <details>
    <summary>Paper</summary>
    Cresswell, Kellen G., and Mikhail G. Dozmorov. "TADCompare: An R Package for Differential and Temporal Analysis of Topologically Associated Domains"  https://doi.org/10.3389/fgene.2020.00158 Frontiers in Genetics 11 (March 10, 2020)
</details>

- [Analysis of the Structural Variability of Topologically Associated Domains as Revealed by Hi-C](https://doi.org/10.1093/nargab/lqz008) - TAD variability among 137 Hi-C samples (including replicates, 69 if not) from 9 studies. HiCrep, Jaccard, TADsim to measure similarity. Variability does not come from genetics. Introduction to TADs. 10-70% of TAD boundaries differ between replicates. 20-80% differ between biological conditions. Much less variation across individuals than across tissue types. Lab -specific source of variation - in situ vs. dilution ligation protocols, restriction enzymes not much. HiCpro to 100kb data, ICE-normalization, [Armatus](#armatus) for TAD calling. Table 1 - all studies and accession numbers.<details>
    <summary>Paper</summary>
    Sauerwald, Natalie, Akshat Singhal, and Carl Kingsford. "Analysis of the Structural Variability of Topologically Associated Domains as Revealed by Hi-C" https://doi.org/10.1093/nargab/lqz008 NAR Genomics and Bioinformatics, 30 September 2019
</details>

- <a name="bpscore">[BPscore](https://github.com/rz6/bp-metric)</a> - metric to compare two TAD segmentations. Formula, methods. More stable to Variation of Information (VI) and Jaccard Index (JI). Python implementation for calculating all three metrics. <details>
    <summary>Paper</summary>
    Zaborowski, Rafał, and Bartek Wilczyński. "BPscore: An Effective Metric for Meaningful Comparisons of Structural Chromosome Segmentations"  https://doi.org/10.1089/cmb.2018.0162 Journal of Computational Biology 26, no. 4 (April 2019)
</details>

- [Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types](https://doi.org/10.1093/bioinformatics/bty265) - Analysis of TAD similarity using variation of information (VI) metric as a local distance measure. Defining structurally similar and variable regions. Comparison with previous studies of genomic similarity. Cancer-normal comparison - regions containing pan-cancer genes are structurally conserved in normal-normal pairs, not in cancer-cancer. [Kingsford-Group/localtadsim](#localtadism). 23 human Hi-C datasets, [Hi-C Pro](#hic-pro) processed into 100kb matrices, [Armatus](#armatus) to call TADs. <details>
    <summary>Paper</summary>
    Sauerwald, Natalie, and Carl Kingsford. "Quantifying the Similarity of Topological Domains across Normal and Cancer Human Cell Types" https://doi.org/10.1093/bioinformatics/bty265 Bioinformatics (Oxford, England), (July 1, 2018)
</details>

- <a name="difftad">[DiffTAD](https://bitbucket.org/rzaborowski/differential-analysis)</a> - differential contact frequency in TADs between two conditions. Two - permutation-based comparing observed vs. expected median interactions, and parametric test considering the sign of the differences within TADs. Both tests account for distance stratum. <details>
    <summary>Paper</summary>
    Zaborowski, Rafal, and Bartek Wilczynski. "DiffTAD: Detecting Differential Contact Frequency in Topologically Associating Domains Hi-C Experiments between Conditions"  https://doi.org/10.1101/093625 BioRxiv, January 1, 2016
</details>

## Prediction of 3D features

- <a name="corigami">C. Origami</a> - a deep neural network (Sequence encoder+Feature encoder->transformer->decoder, Methods) predict cell type-specific Hi-C matrices using DNA sequence (one-hot encoding), CTCF binding (ChIP-seq), and chromatin accessibility (ATAC-seq) profiles (all critical for best performance). Chromosome-wide predictions by joining predictions across sliding windows (2Mb). Performance evaluation - correlation of insulation scores, over 0.95 Pearson. Outperforms [Akita](#akita). Enables prediction of the effect of genetic perturbations. <details>
    <summary>Paper</summary>
    Tan, Jimin, Javier Rodriguez-Hernaez, Theodore Sakellaropoulos, Francesco Boccalatte, Iannis Aifantis, Jane Skok, David Fenyo, Bo Xia, and Aristotelis Tsirigos. “Cell Type-Specific Prediction of 3D Chromatin Architecture.” Preprint. Genomics, March 7, 2022. https://doi.org/10.1101/2022.03.05.483136.
</details>

- <a name="chinn">[ChINN](https://github.com/caofan/chinn)</a> - chromatin interaction neural network, predicting chromatin interactions from DNA sequence. Trained on CTCF- and RNA PolII-mediated loops, as well as on Hi-C data. Gradient boosting trained on functional annotation, distance, or both as predictors. ChINN - CNN trained on sequence. Convergent CTCF orientation is an important predictor, other motifs complement predictive power. Applied to 6 new chronic lymphocytic leukemia samples, patient-specific interactions, vaildated by Hi-C and 4C. <details>
    <summary>Paper</summary>
    Cao, Fan, Yu Zhang, Yichao Cai, Sambhavi Animesh, Ying Zhang, Semih Can Akincilar, Yan Ping Loh, et al. "Chromatin Interaction Neural Network (ChINN): A Machine Learning-Based Method for Predicting Chromatin Interactions from DNA Sequences"  https://doi.org/10.1186/s13059-021-02453-5 Genome Biology, (December 2021)
</details>

- <a name="dl2021-hi-c">[DL2021_HI-C](https://github.com/koritsky/DL2021_HI-C)</a> - deep-learning approach for increasing Hi-C data resolution by appending additional information about genome sequence. Two algorithms: the image-to-image model (modified after [VEHiCLE](#vehicle)), which enhances Hi-C resolution by itself, and the sequence-to-image model (modified after [Akita](#akita)), which uses additional information about the underlying genome sequence for further resolution improvement. Both models are combined with the simple head model (Figure 4). Details of network architecture, training, testing, validation (5 metrics). Various architecture modifications. [VEHiCLE](#vehicle) by itself performs well. <details>
    <summary>Paper</summary>
    Kriukov, Dmitrii, Mark Zaretckii, Igor Kozlovskii, Mikhail Zybin, Nikita Koritskiy, Mariia Bazarevich, and Ekaterina Khrameeva. "Hi-C Resolution Enhancement with Genome Sequence Data"  https://doi.org/10.1101/2021.10.25.465745 Preprint. Systems Biology, October 26, 2021. 
</details>

- <a name="l-hic-reg">[L-HiC-Reg](https://github.com/Roy-lab/Roadmap_RegulatoryVariation)</a> (Local HiC-Reg) - a Random Forest based regression method to predict high-resolution contact counts in new cell lines, and a network-based framework to identify candidate cell line-specific gene networks targeted by a set of variants from a Genome-wide association study (GWAS). Trained on chromosome-specific 1Mb segments in one cell line to predict in another. 55 cell lines, 10 annotations (CTCF, RAD21, TBP, histone marks, DNAse), imputing missing annotations with [Avocado](https://github.com/jmschrei/avocado). Outperforms [GeneHancer](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=geneHancer) and [JEME](http://yiplab.cse.cuhk.edu.hk/jeme/), Networks for 15 GWASs in all cell lines are [available](https://pages.discovery.wisc.edu/~bbaur/Roadmap_RegulatoryVariation/). <details>
    <summary>Paper</summary>
    Baur, Brittany, Jacob Schreiber, Junha Shin, Shilu Zhang, Yi Zhang, Jun S Song, William Stafford Noble, and Sushmita Roy. "Leveraging Epigenomes and Three-Dimensional Genome Organization for Interpreting Regulatory Variation"  https://doi.org/10.1101/2021.08.29.458098 biorXiv, August 30, 2021
</details>

- <a name="ccip">[CCIP](https://github.com/GaoLabXDU/CCIP)</a> (CTCF-mediated Chromatin Interaction Prediction) - predicting CTCF-mediated convergent and tandem loops with transitivity. Transitivity definition from the network of multiple CTCF-interacting regions, convergent and tandem. Incorporating the GCP (graph connecting probability) score, together with CTCF, RAD21, directional CTCF motif one-hot encoding into random forest. GCP is the most important predictive feature. Compared with [Lollipop](https://github.com/ykai16/Lollipop) and [CTCF-MP](https://github.com/ma-compbio/CTCF-MP). <details>
    <summary>Paper</summary>
    Wang, Weibing, Lin Gao, Yusen Ye, and Yong Gao. "CCIP: Predicting CTCF-Mediated Chromatin Loops with Transitivity"  https://doi.org/10.1093/bioinformatics/btab534 Bioinformatics, 20 July 2021.
</details>

- <a name="compartmap">[compartmap](https://bioconductor.org/packages/compartmap/)</a> - A/B compartment reconstruction from bulk and scRNA-seq/scATAC-seq. Steps: preprocessing and summarizing the single-cell assay data; eBayes shrinkage of summarized features towards a global or local target; computing the shrinkage correlation estimate of summarized features; computing domains via SVD (Methods). Also detects genomic rearrangements. Applied to K562 data. [Bioconductor](https://bioconductor.org/packages/compartmap/), [Tweet](https://twitter.com/biobenkj/status/1395014386545201152?s=20). <details>
    <summary>Paper</summary>
    Johnson, Benjamin K, Jean-Philippe Fortin, Kasper D Hansen, Hui Shen, and Triche Jr. "Compartmap Enables Inference of Higher-Order Chromatin Structure in Individual Cells from ScRNA-Seq and ScATAC-Seq"  https://doi.org/10.1101/2021.05.17.444465 preprint, May 18, 2021
</details>

- <a name="sci">[SCI](https://github.com/TheJacksonLaboratory/sci)</a> - A/B sub-compartment prediction from Hi-C data, graph embedding (adaptation of LINE method, better than HOPE and DeepWalk) and k-means clustering. Five subcompartments, informed by eipgenetic modifications. Compared with compartments predicted by HMM (Rao 2014) and k-means clustering by Yaffe-Tanay 2011. Improved network centrality, clustering metrics, enrichment in relevant epigenomic marks, higher intra-loop vs. inter-loop ratio, coordinated gene expression, TF enrichment. XGBoost/RF (allowing feature selection), neural network to predict subcompartments from epigenomics, methylation, sequence-based data, replication timing. Predictions confirmed by ChIA-PET data. GM12878 and K562 data. Python implementation. [sci-DNN](https://github.com/TheJacksonLaboratory/sci-DNN), [News](https://www.jax.org/news-and-insights/2020/march/unfold-the-3d-genome). <details>
    <summary>Paper</summary>
    Ashoor, Haitham, Xiaowen Chen, Wojciech Rosikiewicz, Jiahui Wang, Albert Cheng, Ping Wang, Yijun Ruan, and Sheng Li. “Graph Embedding and Unsupervised Learning Predict Genomic Sub-Compartments from HiC Chromatin Interaction Data.” Nature Communications 11, no. 1 (December 2020): 1173. https://doi.org/10.1038/s41467-020-14974-x.
</details>

- <a name="hi-chip-ml">[Hi-ChIP-ML](https://github.com/MichalRozenwald/Hi-ChIP-ML)</a> - machine learning models for TAD boundary prediction using epigenomic data (ChIP-seq, Histone-seq signal, 5 and 18 feature sets) in Drosophila. Linear regression with 4 regularization types, gradient boosting, bidirectional (to capture both sides around a boundary) LSTM. Methods describe data construction, loss function, model architectures, machine learning framework. One feature at a time to measure feature importance. Chriz factor is the top predictor, then H3K4me1, H3K27me1. Python, sklearn, keras. [Colab notebook](https://colab.research.google.com/drive/1_VEcL2qaMAmd0XyGWriXnii16Lqrsd3c?usp=sharing). <details>
    <summary>Paper</summary>
    Rozenwald, Michal B, Aleksandra A Galitsyna, Grigory V Sapunov, Ekaterina E Khrameeva, and Mikhail S Gelfand. "A Machine Learning Framework for the Prediction of Chromatin Folding in Drosophila Using Epigenetic Features"  https://doi.org/10.7717/peerj-cs.307 PeerJ Comput Sci. 2020 Nov 30
</details>

- <a name="3dpredictor">[3DPredictor](https://github.com/labdevgen/3DPredictor)</a> - Web tool for prediction of enhnancer-promoter interactions from gene expression, CTCF binding and orientation, distance between interacting loci. Also predict changes in 3D genome organization - genomic rearrangements. Benchmarking of [TargetFinder](https://github.com/shwhalen/targetfinder), its performance is overestimated. Class balance slightly improves performance. Two random chromosomes for validation, the rest - for training (10-90% split). Gradient boosting for prediction, significantly better than Random Forest. Multiple performance metrics. Model trained in one cell type can predict in another. Predict quantitative level of interactions. Other tools - [EP2Vec](https://github.com/wanwenzeng/ep2vec), [CTCF-MP](https://github.com/ma-compbio/CTCF-MP), [HiC-Reg](hic-reg). [Web tool](https://genedev.bionet.nsc.ru/Web_3DPredictor/). <details>
    <summary>Paper</summary>
    Belokopytova, Polina S., Miroslav A. Nuriddinov, Evgeniy A. Mozheiko, Daniil Fishman, and Veniamin Fishman. "Quantitative prediction of enhancer–promoter interactions." Genome research 30, no. 1 (2020): 72-84. https://doi.org/10.1101/gr.249367.119
</details>

- <a name="akita">[Akita](https://github.com/calico/basenji/tree/tf2_hic/manuscripts/akita)</a> - Chromatin interaction prediction from sequence only using CNN. Takes in 1Mb sequence and predicts interactions at 2Kb resolution. Includes Distance between interacting regions included. Allows for understanding the effect of mutations. Tested on Rao Hi-C datasets, Micro-C data. [Basenji](https://github.com/calico/basenji/tree/tf2_hic/) network architecture. 80/10/10 training, validation, test sets. <details>
    <summary>Paper</summary>
    Fudenberg, Geoff, David R Kelley, and Katherine S Pollard. "Predicting 3D Genome Folding from DNA Sequence"  https://doi.org/10.1038/s41592-020-0958-x Nature Methods, 12 October 2020
</details>

- <a name="3d-genome-2.0">[3D-GNOME 2.0](https://3dgnome.cent.uw.edu.pl/)</a> - Web-based tool for modeling 3D structure changes due to structural variants (SVs). Reference data - GM12878 ChIA-PET data (CTCF, RNAPII). Changes are modeled by removing or adding contacts between chromatin interaction anchors depending on SV (deletion, duplication, insertion, inversion). [Bitbucket](https://bitbucket.org/4dnucleome/spatial_chromatin_architecture/). <details>
    <summary>Paper</summary>
    Wlasnowolski, Michal, Michal Sadowski, Tymon Czarnota, Karolina Jodkowska, Przemyslaw Szalaj, Zhonghui Tang, Yijun Ruan, and Dariusz Plewczynski. "3D-GNOME 2.0: A Three-Dimensional Genome Modeling Engine for Predicting Structural Variation-Driven Alterations of Chromatin Spatial Structure in the Human Genome"  https://doi.org/10.1093/nar/gkaa388 Nucleic Acids Research, (July 2, 2020)
</details>

- <a name="hic-reg">[HiC-Reg](https://github.com/Roy-lab/HiC-Reg)</a> - Predicting Hi-C contact counts from one-dimensional regulatory signals (Histone marks, CTCF, RAD21, Tbp, DNAse). Random Forest regression. Feature encoding - distance between two regions, pair-concat, window, multi-cell. Works across chromosomes (some chromosomes are worse than others) and cell lines (Gm12878, K562, Huvec, Hmec, Nhek, can be used to predict interactions on new cell lines). Selection of the most important features using multi-task group LASSO (distance, CTCF, Tbp, H4K20me1, DNAse, others). Predicted contacts correspond well to the original contacts (distance-stratified Pearson correlation), define TADs similar to the originals (Jaccard), define significant contacts ([Fit-Hi-C](#fit-hi-c)) more enriched in CTCF binding. Validated on HBA1 and PAPPA gene promoters. Hi-C normalization doesn't have much effect. <details>
    <summary>Paper</summary>
    Zhang, Shilu, Deborah Chasman, Sara Knaack, and Sushmita Roy. "In Silico Prediction of High-Resolution Hi-C Interaction Matrices"  https://doi.org/10.1038/s41467-019-13423-8 Nature Communications 10, no. 1 (December 2019)
</details>

- <a name="tadboundarydetector">[TADBoundaryDectector](https://doi.org/10.1093/nar/gkz315)</a> - TAD boundary prediction from sequence only using deep learning models. 12 architectures tested, with three convolutional and an LSTM layer performed best. Methods, Implementation in Keras-TensorFlow. Model evaluation using different criteria, 96% accuracy reported. Deep learning outperforms feature-based models, among which Boosted Trees, Random Forest, elastic net logistic regression are the best performers. Data augmentation (aka feature engineering) by randomly shifting TAD boundary regions by some base pairs of length (0-100). Tested on Drosophila data. Github depreciated, code unavailable. <details>
    <summary>Paper</summary>
    Henderson, John, Vi Ly, Shawn Olichwier, Pranik Chainani, Yu Liu, and Benjamin Soibam. "Accurate Prediction of Boundaries of High Resolution Topologically Associated Domains (TADs) in Fruit Flies Using Deep Learning"  https://doi.org/10.1093/nar/gkz315 Nucleic Acids Research, May 3, 2019. 
</details>

- <a name="3depiloop">[3DEpiLoop](https://bitbucket.org/4dnucleome/3depiloop)</a> - prediction of 3D interactions from 1D epigenomic profiles using Random Forest trained on CTCF peaks (histone modifications are the most important predictors and TFBSs). <details>
    <summary>Paper</summary>
    Al Bkhetan, Ziad, and Dariusz Plewczynski. "Three-Dimensional Epigenome Statistical Model: Genome-Wide Chromatin Looping Prediction"  https://doi.org/10.1038/s41598-018-23276-8 Scientific Reports 8, no. 1 (December 2018). 
</details>

- <a name="prismr">[PRISMR](https://pure.mpg.de/rest/items/item_3047455/component/file_3047456/content)</a> - a polymer-based method (strings and binders switch SBS polymer model) to model 3D chromatin folding, to predict enhancer-promoter contacts, and to model the effect of structural variations (deletions, duplications, inversions) on the 3D genome organization. Input - Hi-C data. Infers minimal factors that shape chromatin folding and its equilibrium under the laws of physics, without prior assumptions or tunable parameters. Simulated annealing Monte Carlo optimization procedure that minimizes the distance between the predicted polymer model and the input contact matrix under a Bayesian weighting factor to avoid overfitting. Tested on a EPHA4 locus associated with limb malformations, and more. [Newly generated capture Hi-C of mouse limb buds at embryonic day 11.5, and human skin fibroblasts, GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92294). Public murine CH12-LX cells and human IMR90 cells. Implementation - the [LAMMPS package](https://lammps.sandia.gov/), code not available. <details>
    <summary>Paper</summary>
    Bianco, Simona, Darío G. Lupiáñez, Andrea M. Chiariello, Carlo Annunziatella, Katerina Kraft, Robert Schöpflin, Lars Wittler, et al. "Polymer Physics Predicts the Effects of Structural Variants on Chromatin Architecture"  https://doi.org/10.1038/s41588-018-0098-8 Nature Genetics, (May 2018)
</details>

- [Supervised Learning Method for Predicting Chromatin Boundary Associated Insulator Elements](http://www.worldscientific.com/doi/pdf/10.1142/S0219720014420062) - Predicting TAD boundaries using training data and making new predictions. Bayesian network (BNFinder method), random forest vs. basic k-means clustering, ChromHMM, cdBEST. Using sequence k-mers and ChIP-seq data from modENCODE for prediction - CTCF ChIP-seq performs best. Used [Boruta](https://cran.r-project.org/web/packages/Boruta/index.html) package for feature selection. The Bayesian network performs best. <details>
    <summary>Paper</summary>
    Bednarz, Paweł, and Bartek Wilczyński. "Supervised Learning Method for Predicting Chromatin Boundary Associated Insulator Elements"(http://www.worldscientific.com/doi/pdf/10.1142/S0219720014420062 Journal of Bioinformatics and Computational Biology 12, no. 06 (December 2014)
</details>

## SNP-oriented

- <a name="iregnet3d">[iRegNet3D](http://iregnet3d.yulab.org/index/)</a> - Integrated Regulatory Network 3D (iRegNet3D) is a high-resolution regulatory network comprised of interfaces of all known transcription factor (TF)-TF, TF-DNA interaction interfaces, as well as chromatin-chromatin interactions and topologically associating domain (TAD) information from different cell lines. Goal: SNP interpretation. Input: One or several SNPs, rsIDs, or genomic coordinates. Output: For one or two SNPs, on-screen information of their disease-related info, connection over TF-TF and chromatin interaction networks, and whether they interact in 3D and located within TADs. For multiple SNPs, the same info downloadable as text files. <details>
    <summary>Paper</summary>
    Liang, Siqi, Nathaniel D. Tippens, Yaoda Zhou, Matthew Mort, Peter D. Stenson, David N. Cooper, and Haiyuan Yu. "IRegNet3D: Three-Dimensional Integrated Regulatory Network for the Genomic Analysis of Coding and Non-Coding Disease Mutations"  https://doi.org/10.1186/s13059-016-1138-2 Genome Biology 18, no. 1 (December 2017)
</details>

- <a name="hugin">[HUGIn](http://yunliweb.its.unc.edu/HUGIn/)</a> - tissue-specific Hi-C linear display of anchor position and around. Overlay gene expression and epigenomic data. Association of SNPs with genes based on Hi-C interactions. Tissue-specific. <details>
    <summary>Paper</summary>
    Martin, Joshua S, Zheng Xu, Alex P Reiner, Karen L Mohlke, Patrick Sullivan, Bing Ren, Ming Hu, and Yun Li. "HUGIn: Hi-C Unifying Genomic Interrogator"  https://doi.org/10.1093/bioinformatics/btx359 Edited by Inanc Birol. Bioinformatics 33, no. 23 (December 1, 2017)
</details>

- <a name="3dsnp">[3DSNP](http://cbportal.org/3dsnp/)</a> -  3DSNP database integrating SNP epigenomic annotations with chromatin loops. Linear closest gene, 3D interacting gene, eQTL, 3D interacting SNP, chromatin states, TFBSs, conservation. For individual SNPs. <details>
    <summary>Paper</summary>
    Lu, Yiming, Cheng Quan, Hebing Chen, Xiaochen Bo, and Chenggang Zhang. "3DSNP: A Database for Linking Human Noncoding SNPs to Their Three-Dimensional Interacting Genes"  https://doi.org/10.1093/nar/gkw1022 Nucleic Acids Research 45, no. D1 (January 4, 2017)
</details>

## CNV and Structural variant detection

- [EagleC](https://github.com/XiaoTaoWang/EagleC) - structural variant (SVs) detection at high resolution from Hi-C, HiChIP, ChIA-PET, single-cell Hi-C. Unique visual patterns of deletions, duplications, inversions in various read orientations, detected by CNN coupled with ensemble learning (50 models). Detects both intra- and inter-chromosomal SVs. Outperforms Hi-C breakfinder, HiCtrans, HiNT-TL in precision and recall, tested on three breast cancer cells/data (Hi-C and WGS), 26 additional cancer cell lines. Performs well at low sequencing depth, detects fusion genes. Captures SVs and fusion genes missed by WGS and Nanopore sequencing. SV breakpoints occur near TAD boundaries, preferentially form between A-A compartments. ICE-normalized matrices with distance effect removed, Gaussian filter, min-max scaling. Construction of ground truth calls, addressing class imbalance. Delly, smoowe for WGS analysis, sniffles, [Picky](https://github.com/TheJacksonLaboratory/Picky), svim for nanopore analysis, [Arriba](https://arriba.readthedocs.io/en/latest/workflow) for gene fusion detection, [cooltools](https://github.com/open2c/cooltools) for A/B compartment analysis, phased by gene density, HiTAD for TAD detection. [Supplementary material](https://www.science.org/doi/10.1126/sciadv.abn9215): Table S1 - A list of high-confidence SVs for training EagleC models. Table S2 - cancer Hi-C datasets. Table S3 - HiChIP/ChIA-PET datasets. Table S4 - Breakpoint coordinates. Table S5 - Hi-C datasets in normal cell lines or tissues. The list of cancer-related genes was obtained from the [Bushman Lab](http://bushmanlab.org/assets/doc/allOnco_May2018.tsv). [Code and documentation at Zenodo](https://doi.org/10.5281/zenodo.6482060). <details>
  <summary>Paper</summary>
  Wang, Xiaotao, Yu Luan, and Feng Yue. “EagleC: A Deep-Learning Framework for Detecting a Full Range of Structural Variations from Bulk and Single-Cell Contact Maps.” Science Advances 8, no. 24 (June 17, 2022): eabn9215. https://doi.org/10.1126/sciadv.abn9215.
</details>

- <a name="hicpipe">[hicpipe](https://github.com/ChenFengling/HiCpipe)</a> - Alternatively, [HiCnorm](#hicnorm). Normalization preserves CNVs in Hi-C data. <details>
    <summary>Paper</summary>
    Yang, L., Chen, F., Zhu, H., Chen, Y., Dong, B., Shi, M., ... & Zhang, M. Q. (2020). 3D Genome Analysis Identifies Enhancer Hijacking Mechanism for High-Risk Factors in Human T-Lineage Acute Lymphoblastic Leukemia. bioRxiv. 
</details>

- <a name="hint">[HiNT](https://github.com/parklab/HiNT)</a> - CNV and translocation detection from \~10-20% ambiguous chimeric reads in Hi-C data. Three tools: HiNT-Pre - preprocessing of Hi-C data; HiNT-CNV and HiNT-TL - CNV and translocation detection, respectively (accept [HiC-Pro](#hic-pro) output). Tested on K562 (cancer) and Gm12878 (normal) data. Removal of known biases using a GAM with Poisson function. Outperforms [Delly](https://github.com/dellytools/delly), [Meerkat](http://compbio.med.harvard.edu/Meerkat/), [hic_breakfinder](#hic-breakfinder), [HiCtrans](#hictrans). Relatively little overlap with CNVs from WGS (BIC-seq2). Gold-standard - FISH data from Dixon et al., “Integrative Detection and Analysis of Structural Variation in Cancer Genomes.” <details>
    <summary>Paper</summary>
    Wang, Su, Soohyun Lee, Chong Chu, Dhawal Jain, Geoff Nelson, Jennifer M. Walsh, Burak H. Alver, and Peter J. Park. "HiNT: A Computational Method for Detecting Copy Number Variations and Translocations from Hi-C Data"  https://doi.org/10.1101/657080 Preprint. Bioinformatics, June 3, 2019
</details>

- <a name="hic-breakfinder">[hic_breakfinder](https://github.com/dixonlab/hic_breakfinder)</a> - Detection of structural variants (SV) by integrating optical mapping, Hi-C, and WGS. Custom pipeline using [LUMPY](https://github.com/arq5x/lumpy-sv), [Delly](https://github.com/dellytools/delly), [Control-FREEC](https://github.com/BoevaLab/FREEC) software. New Hi-C data on 14 cancer cell lines and 21 previously published datasets. Integration of the detected SVs with genomic annotations, including replication timing. Supplementary data with SVs resolved by individual methods and integrative approaches. <details>
    <summary>Paper</summary>
    Dixon, Jesse R., Jie Xu, Vishnu Dileep, Ye Zhan, Fan Song, Victoria T. Le, Galip Gürkan Yardımcı, et al. "Integrative Detection and Analysis of Structural Variation in Cancer Genomes"  https://doi.org/10.1038/s41588-018-0195-8 Nature Genetics, September 10, 2018.  
</details>

- <a name="tad-fusion-score">[TAD fusion score](https://github.com/HormozdiariLab/TAD-fusion-score)</a> - quantifying the effect of deletions on Hi-C interactions. Intro about TAD fusion effect on genome structure. TAD fusion score - the expected total number of changes in pairwise genomic interactions as a result of the deletion. TAD fusion events are negatively selected. <details>
    <summary>Paper</summary>
    Huynh, Linh, and Fereydoun Hormozdiari. "Contribution of Structural Variation to Genome Structure: TAD Fusion Discovery and Ranking"  https://doi.org/10.1101/279356 BioRxiv, March 9, 2018. 
</details>

- <a name="hicnv">[HiCnv](https://github.com/ay-lab/HiCnv)</a>, <a name="hictrans">[HiCTrans](https://github.com/ay-lab/HiCtrans)</a> - CNV, translocation calling from Hi-C data. CNV calling using HMM on per-restriction site quantified data and 1D-normalized accounting for low GC-content (<0.2), mappability (<0.5). Translocation calling on inter-chromosomal matrices, binned. <details>
    <summary>Paper</summary>
    Chakraborty, Abhijit, and Ferhat Ay. "Identification of Copy Number Variations and Translocations in Cancer Cells from Hi-C Data"  https://doi.org/10.1093/bioinformatics/btx664 Edited by Christina Curtis. Bioinformatics 34, no. 2 (January 15, 2018)
</details>

- [Arima-SV-Pipeline](https://github.com/ArimaGenomics/Arima-SV-Pipeline) - Docker/Singularity image for running [Structural Variant Detection with Arima Hi-C Technology](https://arimagenomics.com/resources/blog/new-products-structural-variant-detection/). Includes several tools ([HiCUP](#hicup), [hic_breakfinder](#hic-breakfinder), [Juicer](#juicer), [HiCCUPS](#hiccups)). 

## Visualization

- <a name="hicube">[HiCube](https://github.com/wmalab/HiCube)</a> - web application for Hi-C map and 3D structure visualization, along with 1D annotations. Uses [HiGlass](https://github.com/higlass/higlass-server). Input - mcool files, db annotation files for HiGlass (created from bigWig, BED, BEDPE, BedGraph, etc.), 3D genome structure data in the [g3d](https://github.com/vimaec/g3d) format. Runs via Docker. [Supplementary Table 1](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btad154/7085595#supplementary-data) - Visualization tool comparison, HiCube vs. Nucleome Browser, WashU Epigenome Browser, HiGlass, 3D Genome Browser,  TADkit, HiC3d Viewer, HiCPlotter, pyGenomeTracks. <details>
  <summary>Paper</summary>
  Ye, Tiantian, Yangyang Hu, Sydney Pun, and Wenxiu Ma. “HiCube: Interactive Visualization of Multiscale and Multimodal Hi-C and 3D Genome Data.” Edited by Can Alkan. Bioinformatics, March 24, 2023, btad154. https://doi.org/10.1093/bioinformatics/btad154.
</details>

- <a name="sashimipy">[Sashimi.py](https://github.com/ygidtu/sashimi.py)</a> - a Python package for genomics data visualization. Supports scRNA-seq, protein–DNA/RNA interactions, long-reads sequencing data, and Hi-C data (BAM, BED, bigWig, bigBed, GTF, BedGraph, h5, more. Strand-aware). Programming and interactive web interfaces. Output: pdf, png, jpg, svg. [Bioconda](https://anaconda.org/bioconda/sashimi-py), [Pypi](https://pypi.org/project/sashimi.py/). [Documentation](https://sashimi.readthedocs.io/). <details>
  <summary>Paper</summary>
  Zhang, Yiming, Ran Zhou, and Yuan Wang. “Sashimi.Py: A Flexible Toolkit for Combinatorial Analysis of Genomic Data.” Preprint. Bioinformatics, November 3, 2022. https://doi.org/10.1101/2022.11.02.514803.
</details>

- <a name="nucleomebrowser">[4D Nucleome Browser](http://www.nucleome.org)</a> - an integrative and multimodal data navigation platform for 4D Nucleome. Synchronized multi-panel visualization of 1D, 2D, 3D models, and microscopy data. Hosts over 2K genomics and 700 image datasets for human (hg38) and mouse (mm10). Supports bigWig, bigBed, Tabix, hic formats. Functionality to compare two conditions, scatterplot for bigWig comparisons. Accepts external data via [Nucleome Bridge](https://tinyurl.com/nb-bridge). Integrated with [HiGlass](http://vis.nucleome.org/static/apps/higlass/). [Documentation](https://nb-docs.readthedocs.io/en/latest/), [Youtube video tutorials](https://tinyurl.com/nb-video-tutorial). [GitHub](https://github.com/nucleome). Local installation available. [Supplementary Notes](https://www.nature.com/articles/s41592-022-01559-3#Sec1) - detailed description and illustration of functionality. Table S1 - functionality comparison with UCSC, WashU, others. Table S2 - URLs. Table S3 - data integrated in the browser. Six tutorials, including [Single cell analysis tutorial](https://github.com/nucleome/Tutorial-SingeCellHiC). [Gallery](https://gallery.nucleome.org/) of examples. <details>
    <summary>Paper</summary>
    Zhu, X., Zhang, Y., Wang, Y. et al. Nucleome Browser: an integrative and multimodal data navigation platform for 4D Nucleome. Nat Methods 19, 911–913 (21 July 2022). https://doi.org/10.1038/s41592-022-01559-3
</details>

- <a name="hicognition">[HiCognition](https://github.com/gerlichlab/HiCognition)</a> - visual exploration of Hi-C data, testing for association with 1D data. Region set concept (ChIP-seq peaks, genes, TADs, etc.). Widget architecture: 1D average, 2D average, embedding/heterogeneity views, association/enrichment analysis. Input: BED, bigWig, [cooler](#cooler) files. Data and precomputations are stored in MySQL database. Prebuild dataset of common features. Docker container. [Documentation](https://gerlichlab.github.io/hicognition/docs/). <details>
    <summary>Paper</summary>
    Langer, Christoph C. H., Michael Mitter, Roman R. Stocsits, and Daniel W. Gerlich. “HiCognition: A Visual Exploration and Hypothesis Testing Tool for 3D Genomics.” Preprint. Bioinformatics, May 1, 2022. https://doi.org/10.1101/2022.04.30.490134.
</details>

- <a name="plotgardener">[plotgardener](https://github.com/PhanstielLab/plotgardener)</a> - R/Bioconductor package for multi-panel genomic and non-genomic data visualization. 45 functions for plotting and annotating multiple genomic data formats: genome sequences, gene/transcript annotations, chromosome ideograms, signal tracks, GWAS Manhattan plots, genomic ranges (e.g. peaks, reads, contact domains, etc), paired ranges (e.g. paired-end reads, chromatin loops, structural rearrangements, QTLs, etc), and 3D chromatin contact frequencies. Auto-recognizes .bam, .bigWig, .hic formats, converts genomic data into standard R objects (e.g., data.frame, tibble, GInteractions). Includes 29 genomes, more through Bioconductor integration. [Documentation](https://phanstiellab.github.io/plotgardener/), [Tweet 1](https://twitter.com/dphansti/status/1436417383732813829?s=20), [Tweet 2](https://twitter.com/mikelove/status/1437388358754443274?s=20). <details>
    <summary>Paper</summary>
    Kramer, Nicole E, Eric S Davis, Craig D Wenger, Erika M Deoudes, Sarah M Parker, Michael I Love, and Douglas H Phanstiel. "Plotgardener: Cultivating Precise Multi-Panel Figures in R"  https://doi.org/10.1101/2021.09.08.459338 Preprint. Bioinformatics, September 10, 2021. 
</details>

- <a name="coolbox">[CoolBox](https://github.com/GangCaoLab/CoolBox)</a> - a [pyGenomeTracks](#pygenometracks)-based Hi-C data visualization in Jupyter and command line, matplotlib plotting system. Works with .hic and .m/cool input, plus .bed, .pairs etc. [Documentation and examples](https://gangcaolab.github.io/CoolBox/gallery.html). <details>
    <summary>Paper</summary>
    Xu, Weize, Quan Zhong, Da Lin, Guoliang Li, and Gang Cao. "CoolBox: A Flexible Toolkit for Visual Analysis of Genomics Data"  https://doi.org/10.1101/2021.04.15.439923 Preprint. Bioinformatics, April 16, 2021
</details>

- <a name="genova">[GENOVA](https://github.com/dewitlab/GENOVA)</a> - an R package for the most common Hi-C analyses and visualization. Quality control - cis/trans ratio, checking for translocations; A/B compartments - single and differential compartment signal (H3K4me1 for compartment assignment); TADs - insulation score and directionality index; Genomic annotations - ChIP-seq signal, gene information; Distance decay, and differential analysis; Saddle plot; Aggregate region/peak/TAD/loop analysis, and differential analysis; Aggregated signal (tornado) plots; Intra-inter-TAD interaction differences; PE-SCAn (paired-end Spatial Chromatin Analysis) and C-SCAn enhancer-promoter aggregation. Data for statistical tests can be extracted from the discovery objects. Applied to Hi-C data from Hap1 cells, WT and deltaWAPL (published) and knockdown cohesin-SA1/SA2 conditions ([HiC-pro](#hic-pro), hg19, [HiCCUPS](#hiccups)). Input - [HiC-pro](#hic-pro), [Juicer](#juicer), .cool files. Storage - compressed sparse triangle format and the user-added metadata. All tools return the "discovery object". [Vignette](https://github.com/robinweide/GENOVA/blob/master/vignettes/GENOVA_vignette.pdf). <details>
    <summary>Paper</summary>
    Weide, Robin H. van der, Teun van den Brand, Judith H.I. Haarhuis, Hans Teunissen, Benjamin D. Rowland, and Elzo de Wit. "Hi-C Analyses with GENOVA: A Case Study with Cohesin Variants"  https://doi.org/10.1101/2021.01.22.427620 Preprint. Genomics, January 24, 2021.
</details>

- <a name="circhic">[circHiC](https://github.com/TrEE-TIMC/circHiC)</a> - circular representation of Hi-C data, projection into polar coordinates. Linear relationship between the radius of a circle and the corresponding genomic distance. Predominantly for bacterial chromosomes. Python implementation. [Documentation](https://tree-timc.github.io/circhic/). <details>
    <summary>Paper</summary>
    Junier, Ivan, and Nelle Varoquaux. "CircHiC: Circular Visualization of Hi-C Data and Integration of Genomic Data"  https://doi.org/10.1101/2020.08.13.249110 Preprint. Bioinformatics, August 14, 2020. 
</details>

- <a name="hicbricks">[HiCBricks](https://bioconductor.org/packages/HiCBricks/)</a> - data format and visualization package. hdf5-based data storage format to handle large Hi-C matrices. Visualization of one or two Hi-C matrices, adding annotations. <details>
    <summary>Paper</summary>
    Pal, Koustav, Ilario Tagliaferri, Carmen M Livi, and Francesco Ferrari. "HiCBricks: Building Blocks for Efficient Handling of Large Hi-C Datasets"  https://doi.org/10.1093/bioinformatics/btz808 Edited by Inanc Birol. Bioinformatics, November 7, 2019
</details>

- <a name="hiceekr">[HiCeekR](https://github.com/lucidif/HiCeekR)</a> - Shiny app and GUI for Hi-C data analysis and interpretation. Input - aligned BAM file, with marked duplicates, restriction enzyme cutting sites (HRF5), genome in FASTQ, optionally ChIP-seq BAM, or RNA-seq gene expression (TSV). The workflow includes filtering (PCR artifacts, self-circle, dangling end fragments, using diffHiC) with diagnostic plots, binning interaction matrices in BED (coordinates) and TSV (counts) formats, normalization (ICE, WavSiS, using chromoR), calling A/B compartments (PCA, using [HiTC](#hitc)), TADs (directionality index, [TopDom](#topdom), [HiCseg](#hicseg)), gene expression/epigenomic integration, network analysis and enrichment in GO, KEGG, other databases (using gProfileR). Visualization of zoomable heatmaps, networks (ggplot2, plotly, heatmaply, networkD3). Starts with creating configuration file. Compared with [GITAR](#gitar), [HiC-Pro](#hic-pro), [HiC-bench](#hic-bench), [HiCdat](#hicdat), [HiCexplorer](#hicexplorer), not [Juicer](#juicer) or [HiGlass](#higlass). Illustrated using Rao2014 Gm12878 data. 32Gb RAM (minimal 16Gb) is sufficient, preprocessing of BAM files (Hi-C or ChIP-seq) is the longest. <details>
    <summary>Paper</summary>
    Di Filippo, Lucio, Dario Righelli, Miriam Gagliardi, Maria Rosaria Matarazzo, and Claudia Angelini. "HiCeekR: A Novel Shiny App for Hi-C Data Analysis"  https://doi.org/10.3389/fgene.2019.01079 Frontiers in Genetics 10 (November 4, 2019): 1079. 
</details>

- [Hi-C data visualization review](https://dev.biologists.org/highwire/markup/1255595/expansion?width=1000&height=500&iframe=true&postprocessors=highwire_tables%2Chighwire_reclass%2Chighwire_figures%2Chighwire_math%2Chighwire_inline_linked_media%2Chighwire_embed) - Good introduction into the 3D genome organization, 115 key references. Table 2. Hi-C visualization tools. <details>
    <summary>Paper</summary>
    Ing-Simmons, Elizabeth, and Juan M. Vaquerizas. "Visualising Three-Dimensional Genome Organisation in Two Dimensions"  https://doi.org/10.1242/dev.177162 Development 146, no. 19 (October 1, 2019)
</details>

- <a name="dnarchitect">[DNARchitect](https://github.com/alosdiallo/DNA_Rchitect)</a> - a Shiny App for visualizing genomic data (HiC, mRNA, ChIP, ATAC, etc.) in bed, bedgraph, and bedpe formats. Wraps [Sushi R package](https://bioconductor.org/packages/release/bioc/html/Sushi.html). [Web-version](http://shiny.immgen.org/DNARchitect/). <details>
    <summary>Paper</summary>
    Ramirez, R N, K Bedirian, S M Gray, and A Diallo. "DNA Rchitect: An R Based Visualizer for Network Analysis of Chromatin Interaction Data"  https://doi.org/10.1093/bioinformatics/btz608 Edited by John Hancock. Bioinformatics, August 2, 2019
</details>

- <a name="higlass">[HiGlass](https://github.com/higlass/higlass)</a> - visualization server for Google maps-style navigation of Hi-C maps. Overlay genes, epigenomic tracks. "Composable linked views" synchronizing exploration of multiple Hi-C datasets linked by location/zoom. [Documentation](https://docs.higlass.io/), [links and resources](https://higlass.io/about). <details>
    <summary>Paper</summary>
    Kerpedjiev, Peter, Nezar Abdennur, Fritz Lekschas, Chuck McCallum, Kasper Dinkla, Hendrik Strobelt, Jacob M. Luber, et al. "HiGlass: Web-Based Visual Exploration and Analysis of Genome Interaction Maps"  https://doi.org/10.1186/s13059-018-1486-1 Genome Biology, (December 2018).
</details>

- <a name="3d-genome-browser">[3D Genome Browser](http://promoter.bx.psu.edu/hi-c/)</a> - visualizing existing Hi-C and other chromatin conformation capture data. Alongside with genomic and epigenomic data. Own data can be submitted in BUTLR format. <details>
    <summary>Paper</summary>
    Wang, Yanli, Fan Song, Bo Zhang, Lijun Zhang, Jie Xu, Da Kuang, Daofeng Li, et al. "The 3D Genome Browser: A Web-Based Browser for Visualizing 3D Genome Organization and Long-Range Chromatin Interactions"  https://doi.org/10.1186/s13059-018-1519-9 Genome Biology 19, no. 1 (December 2018)
</details>

- <a name="hipiler">[HiPiler](https://github.com/flekschas/hipiler)</a> - exploration and comparison of loops and domains as snippets-heatmaps of data. [Documentation](http://hipiler.higlass.io/docs). <details>
    <summary>Paper</summary>
    Lekschas, Fritz, Benjamin Bach, Peter Kerpedjiev, Nils Gehlenborg, and Hanspeter Pfister. "HiPiler: Visual Exploration of Large Genome Interaction Matrices with Interactive Small Multiples"  https://doi.org/10.1109/TVCG.2017.2745978 IEEE Transactions on Visualization and Computer Graphics 24, no. 1 (January 2018). TechBlog: [HiPiler simplifies chromatin structure analysis](http://blogs.nature.com/naturejobs/2017/09/11/techblog-hipiler-simplifies-chromatin-structure-analysis/)
</details>

- <a name="nat">[NAT](https://github.com/laseaman/4D_Nucleome_Analysis_Toolbox)</a> - the 4D Nucleome Analysis Toolbox, for Hi-C data (text, cool format) normalization (ICE, Toeplitz, CNV-Toeplitz), TAD calling (Directionality index, [Armatus](#armatus), custom), karyotype abnormalities visualization on inter-chromosomal matrices, time-course visualization. Matlab. <details>
    <summary>Paper</summary>
    Seaman, Laura, and Indika Rajapakse. "4D Nucleome Analysis Toolbox: Analysis of Hi-C Data with Abnormal Karyotype and Time Series Capabilities"  https://doi.org/10.1093/bioinformatics/btx484 Bioinformatics (Oxford, England) 34, no. 1 (01 2018)
</details>

- <a name="gcmapexplorer">[gcMapExplorer](https://rjdkmr.github.io/gcMapExplorer/)</a> - Genome contact map explorer. Analyze and compare Hi-C contact maps, plot other data types alongside. Normalization: KR, IC, and their distance-centric normalization MCFS (median contact frequency). gcmap and HDF5 file format, ccmap is in numpy memmap forman for fast data access, utilities for conversion of coo (sparse 3-columns), pairCoo, homer, hic, bed, bigWig to gcmap/h5 formats. Command line, GUI, API. Python3. [Documentation](https://gcmapexplorer.readthedocs.io). <details>
  <summary>Paper</summary>
  Kumar, Rajendra, Haitham Sobhy, Per Stenberg, and Ludvig Lizana. “Genome Contact Map Explorer: A Platform for the Comparison, Interactive Visualization and Analysis of Genome Contact Maps.” Nucleic Acids Research 45, no. 17 (September 29, 2017): e152–e152. https://doi.org/10.1093/nar/gkx644.
</details>

- <a name="hic-3dviewer">[HiC-3DViewer](https://bitbucket.org/nadhir/hic3dviewer/src/master/)<a/> - HiC-3DViewer is an interactive web-based tool designed to provide an intuitive environment for investigators to facilitate the 3D exploratory analysis of Hi-C data. It based on Flask and can be run directly or as a docker container. <details>
    <summary>Paper</summary>
    Mohamed Nadhir, Djekidel, Wang, Mengjie, Michael Q. Zhang, Juntao Gao. "HiC-3DViewer: a new tool to visualize Hi-C data in 3D space"  https://doi.org/10.1007/s40484-017-0091-8 Quantitative Biology (2017)
</details>

- <a name="hicplotter">[HiCPlotter](https://github.com/kcakdemir/HiCPlotter)</a> - Hi-C visualization tool, allows for integrating various data tracks. Python implementation. <details>
    <summary>Paper</summary>
    Akdemir, Kadir Caner, and Lynda Chin. "HiCPlotter Integrates Genomic Data with Interaction Matrices"  https://doi.org/10.1186/s13059-015-0767-1 Genome Biology 16 (2015)
</details>

- <a name="nuchart">[NuChart](https://doi.org/10.1371/journal.pone.0075146)</a> - gene-centric network of genes interacting in 3D. Integration of epigenomic features. Statistical network analysis. Code unavailable. <details>
    <summary>Paper</summary>
    Merelli, Ivan, Pietro Liò, and Luciano Milanesi. “NuChart: An R Package to Study Gene Spatial Neighbourhoods with Multi-Omics Annotations.” PloS One 8, no. 9 (2013): e75146. https://doi.org/10.1371/journal.pone.0075146
</details>

- <a name="hitc">[HiTC](https://bioconductor.org/packages/HiTC/)</a> - R package for High Throughput Chromosome Conformation Capture analysis, Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. <details>
    <summary>Paper</summary>
    Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. "HiTC: Exploration of High-Throughput ‘C’ Experiments"  https://doi.org/10.1093/bioinformatics/bts521 Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012)
</details>

- <a name="hicognition">[HiCognition](https://github.com/gerlichlab/HiCognition)</a> - A visual exploration and hypothesis testing tool for 3D genomics. HiCognition is a data exploration tool that allows stream-lined exploration of aggregate genomic data. HiCognition is centered around Hi-C data but also enables integration of Chip-seq and region-based data. Docker installation, uses D3.js and Pixi.js. Similar to [HiGlass](#higlass), [HiPiler](#hipiler). [Documentation](https://gerlichlab.github.io/hicognition/docs/)

- [MuGVRE](https://www.multiscalegenomics.eu/MuGVRE/) - The MuG Virtual Research Environment supports the expanding 3D/4D genomics community by developing tools to integrate and visualize genomics data from sequence to 3D/4D chromatin dynamics data

- <a name="pygenometracks">[pyGenomeTracks](https://github.com/deeptools/pyGenomeTracks)</a> - python module to plot beautiful and highly customizable genome browser tracks.

- <a name="tadkit">[TADKit](https://github.com/3DGenomes/TADkit)</a> - 3D Genome Browser. [Main web site](http://sgt.cnag.cat/3dg/tadkit/).



## De novo genome scaffolding

- Benchmarking of three genome assemblers from Hi-C data (3d-dna, SALSA2, YaHS) on a de novo assembly of Arabidopsis Thaliana. Oxford Nanopore read processing for making the draft assembly. Hi-C scaffolding, overview of each scaffolder, processing with [Arima pipeline](https://github.com/ArimaGenomics/mapping_pipeline), QUAST and BUSCO metrics for benchmarking. YaHS performs best, easiest to install and use. <details>
  <summary>Paper</summary>
  Obinu, Lia, Urmi Trivedi, and Andrea Porceddu. “Benchmarking of Hi-C Tools for Scaffolding de Novo Genome Assemblies.” Preprint. Genomics, May 18, 2023. https://doi.org/10.1101/2023.05.16.540917.
</details>

- Benchmark of five Hi-C scaffolders ([Lachesis](#lachesis), [HiRise](https://github.com/DovetailGenomics/HiRise_July2015_GR), [3D-DNA](#theaidenlab3ddna), [SALSA](https://github.com/marbl/SALSA), [ALLHiC](https://github.com/tanghaibao/allhic)). Accuracy measured by matching scaffolds with the assembly contigs. On average, HiRise and Lachesis performed the best, with HiRise and Salsa working best on less fragmented assemblies, and HiRise, Lacheis, or AllHiC being better choices for more fragmented assemblies. Details and problems with some software. [Docker images](https://hub.docker.com/u/aakashsur) for individual tools. <details>
  <summary>Paper</summary>
  Sur, Aakash, William Stafford Noble, and Peter J. Myler. "A benchmark of Hi-C scaffolders using reference genomes and de novo assemblies." bioRxiv (April 20, 2022). https://doi.org/10.1101/2022.04.20.488415
</details>

- <a name="endhic">[EndHiC](https://github.com/fanagislab/EndHiC)</a> - chromosome scaffolding using Hi-C links from contig ends. Requires contigs from [PacBio's HiFiasm technology](https://github.com/chhylp123/hifiasm) (assembled by [HiCanu](https://github.com/marbl/canu)), and Hi-C data processed by [HiC-Pro](#hic-pro). Applied to human, rice, Arabidopsis, achieves higher accuracy than [Lachesis](#lachesis), [ALLHiC](https://github.com/tanghaibao/allhic), [3D-DNA](#theaidenlab3ddna). Perl scripts. <details>
    <summary>Paper</summary>
    Wang, Sen, Hengchao Wang, Fan Jiang, Anqi Wang, Hangwei Liu, Hanbo Zhao, Boyuan Yang, Dong Xu, Yan Zhang, and Wei Fan. "EndHiC: Assemble Large Contigs into Chromosomal-Level Scaffolds Using the Hi-C Links from Contig Ends](https://arxiv.org/abs/2111.15411 ArXiv, 30 Nov 2021
</details>

- <a name="instagraal">[instaGRAAL](https://github.com/koszullab/instaGRAAL)</a> - reimplementation of [GRAAL genome assembler](#graal) (chromosome level) for large genomes. Similar MCMC approach, implemented on NVIDIA GPU. Tested, among others, on segments of the human genome. <details>
    <summary>Paper</summary>
    Baudry, Lyam, Nadège Guiglielmoni, Hervé Marie-Nelly, Alexandre Cormier, Martial Marbouty, Komlan Avia, Yann Loe Mie, et al. "InstaGRAAL: Chromosome-Level Quality Scaffolding of Genomes Using a Proximity Ligation-Based Scaffolder"  https://doi.org/10.1186/s13059-020-02041-z Genome Biology 21, no. 1 (December 2020)
</details>

- <a name="bin3c">[bin3C](https://github.com/cerebis/bin3C)</a> - resolving metagenome-assembled genomes from Hi-C data. Metagenomic assembly using [SPAdes](http://cab.spbu.ru/software/spades/). Tested using simulated ([Sim3C](https://github.com/cerebis/sim3C) and [MetaART](https://github.com/cerebis/meta-sweeper)) and real-life data. Performance metrics: adjusted mutual information, weighted Bcubed. Contact matrix where bins are contigs. Infomap method for clustering the whole-contig graph. Compared with [ProxiMeta (Phase Genomics)](https://phasegenomics.com/products/proximeta/). <details>
    <summary>Paper</summary>
    DeMaere, Matthew Z., and Aaron E. Darling. "Bin3C: Exploiting Hi-C Sequencing Data to Accurately Resolve Metagenome-Assembled Genomes"  https://doi.org/10.1186/s13059-019-1643-1 Genome Biology 20, no. 1 (December 2019)
</details>

- <a name="hicassembler">[HiCAssembler](https://github.com/maxplanck-ie/HiCAssembler)</a> - Hi-C scaffolding tool combining assembly using Hi-C data with scaffolds from regular sequencing (short or long sequencing). Uses strategies from [Lachesis](#lachesis) and  [3D-DNA](#theaidenlab3ddna). Visual adjustment of scaffolding errors. Automatic and manual misassembly correction. <details>
    <summary>Paper</summary>
    Renschler, Gina, Gautier Richard, Claudia Isabelle Keller Valsecchi, Sarah Toscano, Laura Arrigoni, Fidel Ramirez, and Asifa Akhtar. "Hi-C Guided Assemblies Reveal Conserved Regulatory Topologies on X and Autosomes despite Extensive Genome Shuffling"  https://doi.org/10.1101/580969 BioRxiv, March 18, 2019. 
</details>

- <a name="theaidenlab3ddna">[3D-DNA](https://github.com/theaidenlab/3D-DNA)</a> Hi-C genome assembler and its application/validation. Methods are in the supplemental. [DNA Zoo](https://www.dnazoo.org/methods) - genome assemblies using Hi-C, methods, papers. <details>
    <summary>Paper</summary>
    Dudchenko, Olga, Sanjit S. Batra, Arina D. Omer, Sarah K. Nyquist, Marie Hoeger, Neva C. Durand, Muhammad S. Shamim, et al. "De Novo Assembly of the Aedes Aegypti Genome Using Hi-C Yields Chromosome-Length Scaffolds"  https://doi.org/10.1126/science.aal3327 Science (New York, N.Y.) 356, no. 6333 (07 2017)
</details> <details>
    <summary>Genome Assembly Cookbook</summary>
    https://github.com/theaidenlab/Genome-Assembly-Cookbook) - Genome assembly from Hi-C data, pipeline and instructions from the Aiden Lab. https://aidenlab.org/assembly/manual_180322.pdf
</details>

- <a name="graal">[GRAAL](https://github.com/koszullab/GRAAL)</a> - Genome (Re)Assembly Assessing Likelihood - genome assembly from Hi-C data. Gaps in genome assembly that can be filled by scaffolding. Superior to [Lachesis](#lachesis) and [dnaTri](#dnatri), which are sensitive to duplications, clustering they use to initially arrange the scaffolds, parameters, unknown reliability. A Bayesian approach, prior assumptions are that cis-contact probabilities follow a power-law decay and that counts in the interaction matrix are Poisson. Multiple genomic structures tested using MCMC (Multiple-Try Metropolis algorithm) to maximize the likelihood of data given a genomic structure. <details>
    <summary>Paper</summary>
    Marie-Nelly, Hervé, Martial Marbouty, Axel Cournac, Jean-François Flot, Gianni Liti, Dante Poggi Parodi, Sylvie Syan, et al. "High-Quality Genome (Re)Assembly Using Chromosomal Contact Data"  https://doi.org/10.1038/ncomms6695 Nature Communications 5 (December 17, 2014)
</details>

- <a name="dnatri">[dnaTri](https://github.com/NoamKaplan/dna-triangulation)</a> - genome scaffolding via probabilistic modeling using two constraints of Hi-C data - distance-dependent decay and cis-trans ratio. Using known chromosome scaffolds and de novo assembly. Naive Bayes classifier to distinguish chromosome-specific vs. on different chromosomes contigs. Average linkage clustering to assemble contigs into 23 groups of chromosomes. Completed 65 previously unplaced contigs. [Data](http://my5c.umassmed.edu/triangulation/). <details>
    <summary>Paper</summary>
    Kaplan, Noam, and Job Dekker. "High-Throughput Genome Scaffolding from in Vivo DNA Interaction Frequency"  https://doi.org/10.1038/nbt.2768 Nature Biotechnology 31, no. 12 (December 2013)
</details>

- <a name="lachesis">[Lachesis](https://github.com/shendurelab/LACHESIS)</a> - a three-step genome scaffolding tool: 1) graph clustering of scaffolds to chromosome groups, 2) ordering clustered scaffolds (minimum spanning tree, reassembling longest-to-shortest branches), 3) assigning orientation (exact position and the decay of interactions). Duplications and repeat regions may be incorrectly ordered/oriented. Tested on a normal human, mouse, drosophila genomes, and on the HeLa cancer genome. <details>
    <summary>Paper</summary>
    Burton, Joshua N., Andrew Adey, Rupali P. Patwardhan, Ruolan Qiu, Jacob O. Kitzman, and Jay Shendure. "Chromosome-Scale Scaffolding of de Novo Genome Assemblies Based on Chromatin Interactions"  https://doi.org/10.1038/nbt.2727 Nature Biotechnology 31, no. 12 (December 2013)
</details>

- <a name="yahs">[YaHS](https://github.com/c-zhou/yahs)</a> - yet another Hi-C scaffolding tool. It relies on a new algothrim for contig joining detection which considers the topological distribution of Hi-C signals aiming to distingush real interaction signals from mapping nosies. Implemented in C, very fast.


## 3D modeling

- [The Nucleome Data Bank](https://ndb.rice.edu/) (NDB) - A web platform to simulate and browse the three-dimensional architecture of genomes. Data download from different studies, and computationally simulated. Interactive visualization of these data. Molecular Dynamics simulation using the MEGABASE (Maximum Entropy Genomic Annotation from Biomarkers Associated to Structural Ensembles, neural network to predict A/B compartments from 1D epigenomic data) + MiChroM (Minimal Chromatin Model, energy landscape model, [documentation](https://open-michrom.readthedocs.io/)) computational pipeline. Input: 1D epigenomic data, DNA proximity (3C, HiC, GAM, SPRITE), imaging data. Output: chromosomal folding predictive of HiC, imaging data, GROMACS input files. [.ndb format](https://ndb.rice.edu/ndb-format) to store nucleome structural data, analog of the Protein Data Bank [.pdb](https://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data) format. <details>
  <summary>Paper</summary>
  Contessoto, Vinícius G., Ryan R. Cheng, Arya Hajitaheri, Esteban Dodero-Rojas, Matheus F. Mello, Erez Lieberman-Aiden, Peter G. Wolynes, Michele Di Pierro, and José N. Onuchic. "The Nucleome Data Bank: web-based resources to simulate and analyze the three-dimensional genome." Nucleic Acids Research 49, no. D1 (8 January 2021): D172-D182. https://doi.org/10.1093/nar/gkaa818
</details>

- <a name="pastis-nb">[Pastis-NB](https://github.com/hiclib/pastis)</a> - an extrension of Pastis 3D modeling tool with negative binomial distribution-based modeling. Compared with MDS-based methods ([ShRec3D](#shrec3d), [ChromSDE](https://www.comp.nus.edu.sg/~bioinfo/ChromSDE/), [Pastis-MDS](https://github.com/hiclib/pastis)) and Pastis- PM (Poisson model). More accurate and stable results. [Supplementary material](https://nellev.github.io/pastisnb/). <details>
    <summary>Paper</summary>
    Nelle Varoquaux, William Stafford Noble, Jean-Phillipe Vert. "Inference of genome 3D architecture by modeling overdispersion of Hi-C data"  https://doi.org/10.1101/2021.02.04.429864 bioRxiv. February 05, 2021
</details>

- <a name="taddyn">[TADdyn](https://github.com/3DGenomes/TADbit/tree/TADdyn)</a> - studying time-dependent dynamics of chromatin domains during natural and induced cell processes by simulating smooth 3D transitions of chromosome structure. A part of [TADBit](#tadbit), developed by the Marti-Renom group. Tested on in situ Hi-C time course experiment, reprogramming of murine B cells to pluripotent cells, changes of 21 genomic loci. [Data and video](http://sgt.cnag.cat/3dg/datasets/). <details>
    <summary>Paper</summary>
    Di Stefano, Marco, Ralph Stadhouders, Irene Farabella, David Castillo, François Serra, Thomas Graf, and Marc A. Marti-Renom. "Transcriptional Activation during Cell Reprogramming Correlates with the Formation of 3D Open Chromatin Hubs"  https://doi.org/10.1038/s41467-020-16396-1 Nature Communications 11, no. 1 (December 2020)
</details>

- <a name="stoh-c">[StoH-C](https://github.com/kimmackay/StoHi-C)</a> - 3D genome reconstruction using tSNE. Python scripts for 3D embedding and visualization (plot-ly, matplotlib, Chart Studio). Visually tested on fission yeast genome as compared with MDS-reconstructed genome (wild type, G1-arrested, rad21 mutation, clr4 deletion). <details>
    <summary>Paper</summary>
    MacKay, Kimberly, and Anthony Kusalik. "StoHi-C: Using t-Distributed Stochastic Neighbor Embedding (t-SNE) to Predict 3D Genome Structure from Hi-C Data"  https://doi.org/10.1101/2020.01.28.923615 Preprint. Bioinformatics, January 29, 2020.
</details>

- <a name="hierarchical3dgenome">[Hierarchical3DGenome](https://github.com/BDM-Lab/Hierarchical3DGenome)</a> - high-resolution (5kb) reconstruction of the 3D structure of the genome. Using [LorDG](https://github.com/BDM-Lab/LorDG), first, assemble the 3D model at the level of TADs, then inside individual TADs. Gm12878 cell line, [Arrowhead](#arrowhead) for TAD calling, KR and ICE normalization, benchmarking against [miniMDS](https://github.com/seqcode/miniMDS), five tests including comparison with FISH. <details>
    <summary>Paper</summary>
    Trieu, Tuan, Oluwatosin Oluwadare, and Jianlin Cheng. "Hierarchical Reconstruction of High-Resolution 3D Models of Large Chromosomes"  https://doi.org/10.1038/s41598-019-41369-w Scientific Reports 9, no. 1 (March 21, 2019): 4971. 
</details>

- <a name="csynth">[CSynth](http://csynth.org/)</a> - 3D genome interactive modeling on GPU, and visualization. <details>
    <summary>Paper</summary>
    Todd, Stephen, Peter Todd, Simon J McGowan, James R Hughes, Yasutaka Kakui, Frederic Fol Leymarie, William Latham, and Stephen Taylor. "CSynth: A Dynamic Modelling and Visualisation Tool for 3D Chromatin Structure"  https://doi.org/10.1101/499806 BioRxiv, January 1, 2019
</details>

- [3DMax](https://github.com/BDM-Lab/3DMax) - reconstruction of 3D genome structure. Maximum likelihood objective function. Many simplifying assymptions. Gradient ascent algorithm. Distance Pearson and Spearman correlation coefficients for comparing 3D structures. <details>
  <summary>Paper</summary>
  Oluwadare, Oluwatosin, Yuxiang Zhang, and Jianlin Cheng. “A Maximum Likelihood Algorithm for Reconstructing 3D Structures of Human Chromosomes from Chromosomal Contact Data.” BMC Genomics 19, no. 1 (December 2018). https://doi.org/10.1186/s12864-018-4546-8.
</details>

- <a name="genomeflow">[GenomeFlow](https://github.com/jianlin-cheng/GenomeFlow)</a> - a complete set of tools for Hi-C data alignment, normalization, 2D visualization, 3D genome modeling and visualization. [ClusterTAD](#clustertad) for TAD identification. [LorDG](https://github.com/BDM-Lab/LorDG) and [3DMax](https://github.com/BDM-Lab/3DMax) for 3D genome reconstruction. <details>
    <summary>Paper</summary>
    Trieu, Tuan, Oluwatosin Oluwadare, Julia Wopata, and Jianlin Cheng. "GenomeFlow: A Comprehensive Graphical Tool for Modeling and Analyzing 3D Genome Structure"  https://doi.org/10.1093/bioinformatics/bty802 Bioinformatics (Oxford, England), September 12, 2018. 
</details>

- <a name="shrec3d">[ShRec3D](https://sites.google.com/site/julienmozziconacci/home/softwares)</a> - shortest-path reconstruction in 3D. Genome reconstruction by translating a Hi-C matrix into a distance matrix, then multidimensional scaling. Uses binary contact maps. <details>
    <summary>Paper</summary>
    Lesne, Annick, Julien Riposo, Paul Roger, Axel Cournac, and Julien Mozziconacci. "3D Genome Reconstruction from Chromosomal Contacts"  https://doi.org/10.1038/nmeth.3104 Nature Methods 11, no. 11 (November 2014): 1141–43. 
</details>

- <a name="4dmax">[4DMax](https://github.com/Max-Highsmith/4DMax)</a> - 3D modeling over time, predicts dynamic chromosome conformation using time course Hi-C data. In contrast to [TADdyn](#taddyn), models entire chromosomes, uses gradient descent optimization of a spatial restraint-based maximum-likelihood function. Also in contrast to [TADdyn](#taddyn) that focuses on approx. 2Mb retions, 4DMax models whole chromosomes. Tested on simulate Hi-C progression over 6 time points, and 10-day time course of induced stem cell pluripotency in mice. Preserves and predicts A/B compartments, TADs. Output - video of chromosome dynamics

## Deconvolution

- <a name="thunder">[THUNDER](https://github.com/brycerowland/thundeR)</a> - cell-type deconvolution of Hi-C data. NMF, uses informative interactions within and between chromosomes (top 1000 features by Fano factor), reformatted into matrix form by concatenation. Needs number of cell types k. Tested on *in silico* mixture of cell types. Outperforms [TOAST](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1778-0), [CIBERSORT](https://cibersortx.stanford.edu). R implementation. <details>
    <summary>Paper</summary>
    Rowland, Bryce and Huh, Ruth and Hou, Zoey and Hu, Ming and Shen, Yin and Li, Yun "THUNDER: A Reference-Free Deconvolution Method to Infer Cell Type Proportions from Bulk Hi-C Data"  https://doi.org/10.1101/2020.11.12.379941 bioRxiv, November 12, 2020
</details>

## Haplotype phasing

- <a name="gamibhear">[GAMIBHEAR](https://bitbucket.org/schwarzlab/gamibhear)</a> - an R package for haplotype reconstruction in GAM (Genome Architecture Mapping, nuclear cryosectioned profiles) data. Uses the local proximity of SNVs extended to larger genomic windows using a proximity-scaled graph-based approach. Tested on the F123 mESCs, outperforms [WhatsHap](https://whatshap.readthedocs.io/en/latest/) and [HapCHAT](https://github.com/AlgoLab/HapCHAT), fast. <details>
    <summary>Paper</summary>
    Markowski, Julia, Rieke Kempfer, Alexander Kukalev, Ibai Irastorza-Azcarate, Gesa Loof, Birte Kehr, Ana Pombo, Sven Rahmann, and Roland F Schwarz. “GAMIBHEAR: Whole-Genome Haplotype Reconstruction from Genome Architecture Mapping Data.” Edited by Yann Ponty. Bioinformatics 37, no. 19 (October 11, 2021): 3128–35. https://doi.org/10.1093/bioinformatics/btab238.
</details>

- <a name="haplohic">[HaploHiC](https://github.com/Nobel-Justin/HaploHiC)</a> - Hi-C phasing using SNPs/InDels, placement of unphased reads inferred from the nearby phased reads. <details>
    <summary>Paper</summary>
    Lindsly, Stephen, Wenlong Jia, Haiming Chen, Sijia Liu, Scott Ronquist, Can Chen, Xingzhao Wen, et al. "Functional Organization of the Maternal and Paternal Human 4D Nucleome"  https://doi.org/10.1101/2020.03.15.992164 bioRxiv, June 17, 2021. Supplementary note 2 - algorithm details.
    
- <a name="hichap">[HiCHap](https://pypi.org/project/HiCHap/)</a> - a Python 2.7 package designed to process diploid (and haploid) Hi-C data by using phased SNPs. Propose a novel strategy to correct the systematic biases in diploid Hi-C contact maps, includes VC or ICE normalization. Perform read mapping with read rescue using ligation junction sites, contact map construction based on phased SNPs, whole-genome identification of A/B compartments (PCA_, topological domains (DI) and chromatin loops (HiCCUPS), and allele-specific testing for diploid Hi-C data (permutation, paired t-test, binomial). [GitHub](https://github.com/Prayforhanluo/HiCHap_master). <details>
    <summary>Paper</summary>
    Luo, H., Li, X., Fu, H. et al. HiCHap: a package to correct and analyze the diploid Hi-C data. https://doi.org/10.1186/s12864-020-07165-x BMC Genomics, October 27, 2020
</details>

## Papers

- [3D Genome, from technology to visualization](https://zhonglab.gitbook.io/3dgenome/preface) - a GitBook by Xingzhao Wen and Sheng Zhong covering biological and computational aspects of 3D genomics and RNA-genome interactions.


### Methodological Reviews

- [pipelines_list.csv](pipelines_list.csv) - A list of available pipelines, URLs, from Miura et al., “[Practical Analysis of Hi-C Data](https://link.springer.com/protocol/10.1007%2F978-1-4939-8766-5_16)”

- [pipeline_comparison.csv](pipeline_comparison.csv) - Available analysis options in each pipeline, from Miura et al., “[Practical Analysis of Hi-C Data](https://link.springer.com/protocol/10.1007%2F978-1-4939-8766-5_16)”

- [Table summarizing functionality of Hi-C data analysis tools](https://www.sciencedirect.com/science/article/pii/S1672022918304339?via%3Dihub#t0005), from Calandrelli et al., “GITAR: An Open Source Tool for Analysis and Visualization of Hi-C Data(https://www.sciencedirect.com/science/article/pii/S1672022918304339)”

- [“Analysis Methods for Studying the 3D Architecture of the Genome"](https://doi.org/10.1186/s13059-015-0745-7) - Genome Hi-C technology and methods review. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7#Tab1) - list of tools. Biases, normalization, matrix balancing. Extracting significant contacts, obs/exp ratio, parametric (power-law, neg binomial, double exponential), non-parametric (splines). 3D enrichment. References. TAD identification, directionality index. Outlook, the importance of comparative analysis. <details>
    <summary>Paper</summary>
    Ay, Ferhat, and William S. Noble. “Analysis Methods for Studying the 3D Architecture of the Genome(https://doi.org/10.1186/s13059-015-0745-7).” Genome Biology 16 (September 2, 2015)
</details>

- [“Computational Methods for Assessing Chromatin Hierarchy"](https://doi.org/10.1016/j.csbj.2018.02.003) - Review of higher-order (chromatin conformation capture) and primary order (DNAse, ATAC) technologies and analysis tools. Table 1 - technology summaries. Table 2 - tool summaries. Inter-chromosomal calls using Binarized contact maps. Visualization. Primary order technologies - details and peak calling. <details>
    <summary>Paper</summary>
    Chang, Pearl, Moloya Gohain, Ming-Ren Yen, and Pao-Yang Chen. “Computational Methods for Assessing Chromatin Hierarchy(https://doi.org/10.1016/j.csbj.2018.02.003).” Computational and Structural Biotechnology Journal 16 (2018)
</details>

- [“Computational Methods for Analyzing Genome-Wide Chromosome Conformation Capture Data"](https://doi.org/10.1016/j.copbio.2018.01.023) - 3C-Hi-C tools review, Table 1 lists categorizes main tools, Figure 1 displays all steps in technology and analysis (alignment, resolution, normalization, including accounting for CNVs, A/B compartments, TAD detection, visualization). A concise description of all tools.<details>
    <summary>Paper</summary>
    Nicoletti, Chiara, Mattia Forcato, and Silvio Bicciato. “Computational Methods for Analyzing Genome-Wide Chromosome Conformation Capture Data(https://doi.org/10.1016/j.copbio.2018.01.023).” Current Opinion in Biotechnology 54 (December 2018)
</details>

- [“Hi-C Analysis: From Data Generation to Integration"](https://doi.org/10.1007/s12551-018-0489-1) - Hi-C technology, data, 3D structures, analysis, and tools. Technology improvement and increasing resolution. FASTQ processing steps ("Hi-C data analysis: from FASTQ to interaction maps" section), pipelines, finding minimum resolution, normalization. Downstream analysis: A/B compartment detection, TAD callers, Hierarchical TADs, interaction callers. Data formats (pairix, sparse matrix format, cool, hic, butlr, hdf5, pgl). Hi-C visualization tools. [Table 2 - summary and comparison of all tools](https://link.springer.com/article/10.1007%2Fs12551-018-0489-1#Tab2)<details>
    <summary>Paper</summary>
    Pal, Koustav, Mattia Forcato, and Francesco Ferrari. “Hi-C Analysis: From Data Generation to Integration(https://doi.org/10.1007/s12551-018-0489-1).” Biophysical Reviews, December 20, 2018.
</details>

- [“Software Tools for Visualizing Hi-C Data"](https://doi.org/10.1186/s13059-017-1161-y) - Hi-C technology, data, and visualization review. Suggestions about graph representation.<details>
    <summary>Paper</summary>
    Yardımcı, Galip Gürkan, and William Stafford Noble. “Software Tools for Visualizing Hi-C Data(https://doi.org/10.1186/s13059-017-1161-y).” Genome Biology 18, no. 1 (December 2017).
</details>

- [“Storage, Visualization, and Navigation of 3D Genomics Data"](https://doi.org/10.1016/j.ymeth.2018.05.008) - Review of tools for visualization of 3C-Hi-C data, challenges, analysis (Table 1). Data formats (hic, cool, BUTLR, ccmap). Database to quickly access 3D data. Details of each visualization tool in Section 4. <details>
    <summary>Paper</summary>
    Waldispühl, Jérôme, Eric Zhang, Alexander Butyaev, Elena Nazarova, and Yan Cyr. “Storage, Visualization, and Navigation of 3D Genomics Data(https://doi.org/10.1016/j.ymeth.2018.05.008).” Methods, May 2018
</details>

- [“An Overview of Methods for Reconstructing 3-D Chromosome and Genome Structures from Hi-C Data"](https://doi.org/10.1186/s12575-019-0094-0) - 3D genome reconstruction review. Intro into equilibrium/fractal globule models. Classification of reconstruction methods: distance-, contact-. and probability-based. [Table 1](https://biologicalproceduresonline.biomedcentral.com/articles/10.1186/s12575-019-0094-0#Tab1) summarizes many tools, methods, and references.<details>
    <summary>Paper</summary>
    Oluwadare, Oluwatosin, Max Highsmith, and Jianlin Cheng. “An Overview of Methods for Reconstructing 3-D Chromosome and Genome Structures from Hi-C Data.(https://doi.org/10.1186/s12575-019-0094-0)” Biological Procedures Online 21, no. 1 (December 2019)
</details>


### General Reviews

- [“Getting the Genome in Shape: The Formation of Loops, Domains and Compartments"](https://doi.org/10.1186/s13059-015-0730-1).” - TAD/loop formation review. Convergent CTCF, cohesin, mediator, different scenarios of loop formation. Stability and dynamics of TADs. Rich source of references.<details>
    <summary>Paper</summary>
    Bouwman, Britta A. M., and Wouter de Laat. “Getting the Genome in Shape: The Formation of Loops, Domains and Compartments(https://doi.org/10.1186/s13059-015-0730-1).” Genome Biology 16 (August 10, 2015)
</details>

- [“The Role of 3D Genome Organization in Disease: From Compartments to Single Nucleotides"](https://doi.org/10.1016/j.semcdb.2018.07.005) - 3D genome structure and disease. Evolution of technologies from FISH to variants of chromatin conformation capture. Hierarchical 3D organization, Table 1 summarizes each layer and its involvement in disease. Rearrangement of TADs/loops in cancer and other diseases. Specific examples of the biological importance of TADs, loops as means of distal communication.<details>
    <summary>Paper</summary>
    Chakraborty, Abhijit, and Ferhat Ay. “The Role of 3D Genome Organization in Disease: From Compartments to Single Nucleotides(https://doi.org/10.1016/j.semcdb.2018.07.005).” Seminars in Cell & Developmental Biology 90 (June 2019): 104–13.
</details>

- ["The role of 3D genome organization in development and cell differentiation"](https://www.nature.com/articles/s41580-019-0132-4) - 3D structure of the genome and its changes during gametogenesis, embryonic development, lineage commitment, differentiation. Changes in developmental disorders and diseases. Chromatin compartments and TADs. Chromatin changes during X chromosome inactivation. Promoter-enhancer interactions established during development are accompanied by gene expression changes. Polycomb-mediated interactions may repress developmental genes. References to many studies. <details>
    <summary>Paper</summary>
    Zheng, H., and Xie, W. (2019). "The role of 3D genome organization in development and cell differentiation(https://www.nature.com/articles/s41580-019-0132-4)." Nat. Rev. Mol. Cell Biol.
</details>

- [“The Three-Dimensional Organization of Mammalian Genomes"](https://doi.org/10.1146/annurev-cellbio-100616-060531) - 3D genome structure review. The role of gene promoters, enhancers, and insulators in regulating gene expression. Imaging-based tools, all flavors of chromatin conformation capture technologies. 3D features - chromosome territories, topologically associated domains (TADs), the association of TAD boundaries with replication domains, CTCF binding, transcriptional activity, housekeeping genes, genome reorganization during mitosis. Use of 3D data to annotate noncoding GWAS SNPs. 3D genome structure change in disease.<details>
    <summary>Paper</summary>
    Yu, Miao, and Bing Ren. “The Three-Dimensional Organization of Mammalian Genomes(https://doi.org/10.1146/annurev-cellbio-100616-060531).” Annual Review of Cell and Developmental Biology 33 (06 2017)
</details>

- [“Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation"](http://msb.embopress.org/content/msb/11/12/852.full.pdf) - 3D genome organization parts. Well-written and detailed. References. Technologies: FISH, 3C. 4C, 5C, Hi-C, GCC, TCC, ChIA-PET. Typical resolution - 40bp to 1Mb. LADs - conserved, but some are cell type-specific. Chromosome territories. Cell type-specific. Inter-chromosomal interactions may be important to define cell-specific interactions. A/B compartments identified by PCA. Chromatin loops, marked by CTCF and Cohesin binding, sometimes, with Mediator. Transcription factories. <details>
    <summary>Paper</summary>
    Fraser, J., C. Ferrai, A. M. Chiariello, M. Schueler, T. Rito, G. Laudanno, M. Barbieri, et al. “Hierarchical Folding and Reorganization of Chromosomes Are Linked to Transcriptional Changes in Cellular Differentiation(http://msb.embopress.org/content/msb/11/12/852.full.pdf).” Molecular Systems Biology 11, no. 12 (December 23, 2015)
</details>

- [“Exploring the Three-Dimensional Organization of Genomes: Interpreting Chromatin Interaction Data"](https://doi.org/10.1038/nrg3454) - 3D genome review. Chromosomal territories, transcription factories. Details of each 3C technology. Exponential decay of interaction frequencies. Box 2: A/B compartments (several Mb), TAD definition, size (hundreds of kb). TADs are largely stable, A/B compartments are tissue-specific. Adjacent TADs are not necessarily of opposing signs, may jointly form A/B compartments. Genes co-expression, enhancer-promoters interactions are confined to TADs. 3D modeling. <details>
    <summary>Paper</summary>
    Dekker, Job, Marc A. Marti-Renom, and Leonid A. Mirny. “Exploring the Three-Dimensional Organization of Genomes: Interpreting Chromatin Interaction Data(https://doi.org/10.1038/nrg3454).” Nature Reviews. Genetics 14, no. 6 (June 2013)
</details>

- [“On the Assessment of Statistical Significance of Three-Dimensional Colocalization of Sets of Genomic Elements"](https://doi.org/10.1093/nar/gks012) <details>
    <summary>Paper</summary>
    Witten, Daniela M., and William Stafford Noble. “On the Assessment of Statistical Significance of Three-Dimensional Colocalization of Sets of Genomic Elements(https://doi.org/10.1093/nar/gks012).” Nucleic Acids Research 40, no. 9 (May 2012)
</details>

### Technology

- [4D Nucleome Protocols](https://www.4dnucleome.org/protocols.html) - Collection of genomic technologies currently in use or being developed in the 4DN network - links to wet-lab protocols and papers. [4DN portal blog](https://www.4dnucleome.org/outreach.html)

- [liCHi-C](https://github.com/JavierreLab/liCHiC) - low input capture Hi-C technology (above 50K cells in contrast to 30-50M cells for promoter-capture). Developed for promoter capture but can be adapted to any capture regions. HiCUP pipeline, CHiCAGO to detect significant interactions, HICCUPS, Mustache and HiCExplorer to call loops on Knight-Ruiz normalized matrices at 5kb resolution. Applied to normal and malignant human hematopoietic hierarchy in clinical samples, reconstructs lineages. Benchmarked against Low-C, Hi-C with different restriction enzymes, TagHi-C. Can be used to detect structural variants, breakpoints, link disease variants to genes. [Data](https://ega-archive.org/studies/EGAS00001006305), scripts and processed data on [Zenodo](https://zenodo.org/record/7351026). Additional software: [CHIGP](https://github.com/ollyburren/CHIGP) R package for integrating GWAS summary statistics with CHiCAGO output, poor man's imputation algorithm for GWAS, Blockshifter, COGS. [RegioneReloaded](https://bioconductor.org/packages/regioneReloaded/) - Bioconductor package for comparative permutation analysis of multiple genomic region sets. [HiCaptuRe](https://github.com/LaureTomas/HiCaptuRe/) - R package for Capture Hi-C data management, handles CHiCAGO output. See also the [promoter-capture PCHi-C paper](https://github.com/mdozmorov/HiC_data#phic). [GitHub](https://github.com/JavierreLab/liCHiC/tree/main). <details>
  <summary>Paper</summary>
  Tomás-Daza, Laureano, Llorenç Rovirosa, Paula López-Martí, Andrea Nieto-Aliseda, François Serra, Ainoa Planas-Riverola, Oscar Molina, et al. “Low Input Capture Hi-C (LiCHi-C) Identifies Promoter-Enhancer Interactions at High-Resolution.” Nature Communications 14, no. 1 (January 17, 2023): 268. https://doi.org/10.1038/s41467-023-35911-8.
</details>

- **NOMe-HiC** - 3D, chromatin accessibility, DNA methylation, and genetic variant profiling. Extends [NOME-seq](http://www.genome.org/cgi/doi/10.1101/gr.143008.112) (nucleosome occupancy and methylome sequencing), nucleosome position and methylation measures, uses a GpC methyltransferase, [Methyl-HiC](https://doi.org/10.1038/s41592-019-0502-z), and [Bis-SNP](https://doi.org/10.1186/gb-2012-13-7-r61) that combines DNA methylation and SNP calling from Bisulfite-seq data. Applied to IMR-90, GM12878, high similarity with previous omics measures. Excluded measurements at ambigious CCG sites. Resolution still limited, 20kb currently. Methods with details of analyzing each data modality. Data on [GSE189158](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189158). <details>
  <summary>Paper</summary>
  Fu, Hailu, Haizi Zheng, Xiaoting Chen, Matthew T. Weirauch, Louis J. Muglia, Li Wang, and Yaping Liu. “NOMe-HiC: Joint Profiling of Genetic Variants, DNA Methylation, Chromatin Accessibility, and 3D Genome in the Same DNA Molecule.” Preprint. Genomics, March 30, 2022. https://doi.org/10.1101/2022.03.29.486102.
</details>

- **ChIATAC** technology, combines proximity ligation in ChIA-PET and transposase accessibility in ATAC-seq to identify open chromatin loci and associated chromatin loops. HpyCH4V 4-cutter restriction enzyme, can be combined with AluI. ChIATAC profiles better correlate with RNAPII ChIP-seq, not CTCF, enriched in promoters. Compared with HiCAR, a similar technology, provides more detailed looping structures between active gene promoters and enhancers. Loop calling with HiCCUPS and Mustache, differential TF analysis, TADs calling using HiCExplorer, APA analysis. Applied to Drosophila S2 cells, GM12878, T cell activation by T cell receptor (TCR) and interleukin-2 (IL-2), [GSE194036](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194036). <details>
  <summary>Paper</summary>
  Chai, Haoxi, Harianto Tjong, Peng Li, Wei Liao, Ping Wang, Chee Hong Wong, Chew Yee Ngan, Warren J. Leonard, Chia-Lin Wei, and Yijun Ruan. “ChIATAC Is an Efficient Strategy for Multi-Omics Mapping of 3D Epigenomes from Low-Cell Inputs.” Nature Communications 14, no. 1 (January 13, 2023): 213. https://doi.org/10.1038/s41467-023-35879-5.
</details>

- [pHi-C 3.0 protocol](https://doi.org/10.1101/2020.12.26.424448) - Evaluation of 12 experimental Hi-C protocols - they capture different 3D genome features with different efficiencies. Additional crosslinking with DSG improves signal-to-noise, loop detection, reduced compartment detection. Evaluating 4 restriction enzymes, MNase, DdeI, DpnII, HindIII. 4 cell lines - H1-hESCs, differentiated endoderms, HFF, HeLa (two cell cycle stages). 63 libraries total. All protocols detect cell type-specific differences, A/B compartments, insulation strength. MNase digestion improves loop detection. Anchors for multi-loop interactions can be detected. Double enzyme use improves loop detection. Evaluation of enrichment of CTCF, SMC3, H3K4me3, H3K27ac at loop boundaries. Ultra-deeply sequenced maps using Hi-C, Micro-C, and Hi-C 3.0 protocols (not yet available). [cLIMS](https://github.com/dekkerlab/cLIMS-docker) Hi-C data management system. [Scripts to reproduce results](https://github.com/dekkerlab/matrix_paper). [Tweet](https://twitter.com/job_dekker/status/1343529467277430784?s=20). <details>
    <summary>Paper</summary>
    Oksuz, Betul Akgol, Liyan Yang, Sameer Abraham, Sergey V. Venev, Nils Krietenstein, Krishna Mohan Parsi, Hakan Ozadam, et al. "Systematic Evaluation of Chromosome Conformation Capture Assays"  https://doi.org/10.1101/2020.12.26.424448 Preprint. Genomics, December 27, 2020.
</details>

- [Chromatin conformation capture technologies, from 3C to imaging](https://doi.org/10.1016/j.molcel.2019.12.021) -  3D structures (nucleolus, nuclear speckles, polycomb bodies, chromosome territories, A/B compartments, TADs, loops), their roles in gene expression (sometimes, conflicting), replication timing, DNA repair. Agreement (and disagreement) examples between 3C methods and [FISH](#fish). Description of libation-based (3C - Hi-C) and ligation-free (GAM, SPRITE, DamC) technologies. Multiway interactions, primarily occur within TADs. TAD formation, loop extrusion mechanisms. Association between replication timing and A/B compartments. Effect of mechanical forces on chromosome folding. <details>
    <summary>Paper</summary>
    McCord, Rachel Patton, Noam Kaplan, and Luca Giorgetti. "Chromosome Conformation Capture and Beyond: Toward an Integrative View of Chromosome Structure and Function"  https://doi.org/10.1016/j.molcel.2019.12.021 Molecular Cell, January 2020
</details>

![](/img/3C_technologies.png)

- [Review of Hi-C, Capture-C, and Capture-C technologies, their computational preprocessing](https://doi.org/10.3390/genes10070548) - Experimental protocols, similarities and differences, types of reads (figures), details of alignment, read orientation, elimination of artifacts, quality metrics. A brief overview of preprocessing tools. Example preprocessing of three types of data. Java tool for preprocessing all types of data, [Diachromatic](https://github.com/TheJacksonLaboratory/diachromatic) (Differential Analysis of Chromatin Interactions by Capture), [GOPHER](#gopher) (Generator Of Probes for capture Hi-C Experiments at high Resolution) for genome cutting, probe design. <details>
    <summary>Paper</summary>
    Hansen, Peter, Michael Gargano, Jochen Hecht, Jonas Ibn-Salem, Guy Karlebach, Johannes T. Roehr, and Peter N. Robinson. "Computational Processing and Quality Control of Hi-C, Capture Hi-C and Capture-C Data"  https://doi.org/10.3390/genes10070548 Genes 10, no. 7 (July 18, 2019): 548. 
</details>

- [Chromosome conformation capture technologies, 4C, 5C, Hi-C, ChIP-loop, ChIA-PET](https://doi.org/10.1101/gad.179804.111) - From microscopy observations (constrained movement of genomic loci, LADs, preferential stability of chromosome conformation and its independence from transcription), to technology details (Figure 1). Examples of alpha- and beta-globin locus studies by different technologies, X chromosome inactivation, HOXA-d gene clusters. The future vision of single-cell, single-allele investigation of chromatin interactions. <details>
    <summary>Paper</summary>
    Wit, E. de, and W. de Laat. "A Decade of 3C Technologies: Insights into Nuclear Organization"  https://doi.org/10.1101/gad.179804.111 Genes & Development 26, no. 1 (January 1, 2012)
</details>

- [Review of technologies for studying the 3D structure of the genome](https://doi.org/10.1101/gad.281964.116). From microscopy to 3C techniques revealing CTCF and cohesin as the key proteins for establishing chromatin loops.TADs are unlikely over large distances >>1Mb. Details of 3C, 4C, 5C, Hi-C, ChIP-PET, and other derivatives. A/B compartments and their subdivision. TADs, their conservation, \~35-50% still seem to change. CTCF (directionality of binding important) and cohesin. Diseases and the 3D genome, examples. Key steps in data analysis and interpretation, software, visualization. Hi-C data specifics - chimeric reads, mapping, data representation as fixed or enzyme-sized bins, normalization, detection of A/B compartments and TAD boundaries, significant interactions. Hi-C analysis tools: [HiC-Pro](#hic-pro), [HiCUP](#hicup), [HOMER](#homer), [Juicer](#juicer). Tools for 3D modeling. <details>
    <summary>Paper</summary>
    Denker, Annette, and Wouter de Laat. "The Second Decade of 3C Technologies: Detailed Insights into Nuclear Organization"  https://doi.org/10.1101/gad.281964.116 Genes & Development 30, no. 12 (June 15, 2016)
</details>

- [RedChIP](https://github.com/agalitsyna/RedClib) -  RedClib processing pipeline; technology for identifying RNA-DNA interactions mediated by a particular protein. Combines RNA-DNA proximity ligation and ChIP. EZH2 in human ESCs and CTCF in K562. An input fraction without IP sequenced in parallel. [GEO GSE174474](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&hx0026;acc=GSE174474) - RedChIP replicates (EZH2/CTCF samples and input). See the [RedC paper](https://doi.org/10.1093/nar/gkaa457) for more details. <details>
    <summary>Paper</summary>
    Gavrilov, Alexey A., Rinat I. Sultanov, Mikhail D. Magnitov, Aleksandra A. Galitsyna, Erdem B. Dashinimaev, Erez Lieberman Aiden, and Sergey V. Razin. “RedChIP Identifies Noncoding RNAs Associated with Genomic Sites Occupied by Polycomb and CTCF Proteins.” Proceedings of the National Academy of Sciences 119, no. 1 (January 4, 2022): e2116222119. https://doi.org/10.1073/pnas.2116222119.
</details>

- [RD-SPRITE](https://github.com/GuttmanLab/sprite2.0-pipeline) -  (RNA & DNA SPRITE) technology to map the spatial organization of RNA and DNA (improved efficiency of the RNA-tagging steps). Noncoding RNAs form high-concentration territories throughout the nucleus (ncRNAs involved in rRNA processing organize around rRNA genes, splicing ncRNAs are concentrated near PolII sites, more examples). In general, (noncoding) RNA-chromatin structures are associated with three functional classes: RNA processing, heterochromatin assembly, gene regulation. [GSE GEO151515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151515) - RD-SPRITE data (clusters information) on mouse ESCs, naive and treated with Actinomycin D. See also the [other sprite pipeline](https://github.com/GuttmanLab/sprite-pipeline) <details>
    <summary>Paper</summary>
    Quinodoz, Sofia A., Joanna W. Jachowicz, Prashant Bhat, Noah Ollikainen, Abhik K. Banerjee, Isabel N. Goronzy, Mario R. Blanco, et al. “RNA Promotes the Formation of Spatial Compartments in the Nucleus.” Cell 184, no. 23 (November 2021): 5775-5790.e30. https://doi.org/10.1016/j.cell.2021.10.014.
</details>

- [Hi-CO](https://figshare.com/articles/software/Hi-CO_3D_genome_structure_analysis_with_nucleosome_resolution/13176101/1) - Example data and tutorial of chromatin conformation capture with nucleosome orientation technology. Ligation of DNA entry or exit points at every nucleosome locus in a micrococcal nuclease (MNase)-fragmented genome. Produces \~300M reads. Computational analysis - simulated annealing - molecular dynamics, determines 3D positions and orientation of every nucleosome. Not suitable for single-cell genome analysis, only detects pairwise ligations. Applied to the S. cerevisiae genome. Protocol. <details>
    <summary>Paper</summary>
    - Ohno, Masae, Tadashi Ando, David G. Priest, and Yuichi Taniguchi. "Hi-CO: 3D Genome Structure Analysis with Nucleosome Resolution"  https://doi.org/10.1038/s41596-021-00543-z Nature Protocols, May 28, 2021. 
</details>

- [RADICL-seq](https://github.com/fagostini/RADICL_analysis) - RNA And DNA Interacting Complexes Ligated and sequenced - RNA-chromatin interaction profiling in intact nuclei. Crosslinking, digesting withDNAse I, RNAse H treatment (reduces ribosomal contamination), a bridge adapter specifically ligating RNA and DNA, reverse crosslinking, conversion to double-stranded DNA, cutting to 25-27bp length with EcoP15I enzyme. RNA- and DNA tags mapped, quantified by distance. Improvement over other technologies (MARGI, GRID-seq, SPRITE) - less cells, less spurious interactions, more trans interactions (lncRNAs in particular), more long-distance interactions. Four technology advantages in Discussion. Analysis - accounts for coverage, one-sided binomial test to call significant interactions. Data - mESC, oli-neu (from mOPCs). Visualization as 25kb contact matrix. Several noncoding transcripts interact in trans. Integration with CAGE data, comparison with Hi-C, association with TAD boundaries, A/B compartments, repeat elements. Differential analysis, Cell-type-specific interactions. [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132192). <details>
    <summary>Paper</summary>
    Bonetti, Alessandro, Federico Agostini, Ana Maria Suzuki, Kosuke Hashimoto, Giovanni Pascarella, Juliette Gimenez, Leonie Roos, et al. "RADICL-Seq Identifies General and Cell Type–Specific Principles of Genome-Wide RNA-Chromatin Interactions"  https://doi.org/10.1038/s41467-020-14337-6 Nature Communications, (December 2020)
</details>

- [tagHi-C](https://doi.org/10.1016/j.celrep.2020.108206) - low-input tagmentation-based Hi-C. Applied to mouse hematopoiesis 10 major blood cell types. Changes in compartments and the Rabl configuration defining chromatin condensation. Gene-body-associating domains are a general property of highly-expressed genes. Spatial chromatin loops link GWAS SNPs to candidate blood-phenotype genes. [HiC-Pro](#hic-pro) to [Juicer](#juicer). [GEO GSE142216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142216) - RNA-seq, replicates, [GEO GSE152918](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152918) - tagHi-C data, replicates, combined .hic files. <details>
    <summary>Paper</summary>
    Zhang C, Xu Z, Yang S, Sun G, Jia L, Zheng Z, Gu Q, Tao W, Cheng T, Li C, Cheng H. "tagHi-C Reveals 3D Chromatin Architecture Dynamics during Mouse Hematopoiesis"  https://doi.org/10.1016/j.celrep.2020.108206 Cell Rep. 2020 Sep 29 
</details>

- [scsHi-C](https://github.com/gerlichlab/scshic_analysis) - sister-chromatid-sensitive Hi-C to explore interactions between the sister chromatids. Distinguishing cis from trans sister contacts based on 4-thio-thymidine (4sT) labeling. Paired organization of sister chromatins in interphase and complete separation in mitosis. TADs that exhibit tight pairing are heterochromatin marked by H3K27me3. Chromatids are predominantly linked at TAD boundaries, within TADs - more flexible. Investigation of looping mechanism - NIPBL-depletion, Sororin degradation. Jupyter notebooks for each analysis. <details>
    <summary>Paper</summary>
    Mitter, Michael, Catherina Gasser, Zsuzsanna Takacs, Christoph C. H. Langer, Wen Tang, Gregor Jessberger, Charlie T. Beales, et al. "Sister-Chromatid-Sensitive Hi-C Reveals the Conformation of Replicated Human Chromosomes"  https://doi.org/10.1101/2020.03.10.978148 Preprint. Cell Biology, March 11, 2020. 
</details>

- [4C technology, wet-lab protocol, and data analysis and visualization](https://doi.org/10.1016/j.ymeth.2019.07.014) -  R-based processing pipeline [pipe4C](https://github.com/deLaatLab/pipe4C), configuration parameters. <details>
    <summary>Paper</summary>
    Krijger, Peter H.L., Geert Geeven, Valerio Bianchi, Catharina R.E. Hilvering, and Wouter de Laat. "4C-Seq from Beginning to End: A Detailed Protocol for Sample Preparation and Data Analysis"  https://doi.org/10.1016/j.ymeth.2019.07.014 Methods 170 (January 2020)
</details>

- [SisterC Hi-C technology to test interactions between sister chromatids](http://biorxiv.org/content/early/2020/03/11/2020.03.10.986208.abstract) -  Uses BrdU incorporation in S-phase and single-strand degradation by UV/Hoechst treatment to obtain inter-sister or intra-sister interactions. Findings about the alignment of chromatids, strong at centromeres, looser (\~35kb spaced) interactions along arms, loops up to 50kb. Tested on the S. cerevisiae genome. distiller-nf, pairtools, cooltools for processing. <details>
    <summary>Paper</summary>
    "Detecting Chromatin Interactions along and between Sister Chromatids with SisterC(http://biorxiv.org/content/early/2020/03/11/2020.03.10.986208.abstract bioRxiv 2020.03.10
</details>

- [AQuA-HiChIP](https://doi.org/10.1038/s41596-019-0285-9) - absolute quantification of chromatin architecture HiChIP technology to quantify the absolute differences in protein-anchored chromatin interactions. Wet-lab and bioinformatics protocols. Involves a mixture of human and mouse material in a predefined proportion (e.g., 3:1 should be equal for both samples), used to normalize 3D interactions on exogenous chromatin from an orthogonal species. AQuA normalization by calculating the AQuA factor as the ratio of human to mouse reads and applying this sample-specific scaling factor to its valid contact matrix at a locus of interest. HiC-Pro processing. [GSE120770](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120770) - AQuA-HiChIP data, Genome-wide 3D interaction profiles as mediated by H3K27ac histone marks, in Rhabdomyosarcoma (RMS) cell lines before and after HDAC inhibition. R-based pipeline to generate contact matrices and APA plots with and without AQuA normalization, [GitHub](https://github.com/GryderArt/AQuA-HiChIP). <details>
  <summary>Paper</summary>
  Gryder, Berkley E., Javed Khan, and Benjamin Z. Stanton. “Measurement of Differential Chromatin Interactions with Absolute Quantification of Architecture (AQuA-HiChIP).” Nature Protocols, February 12, 2020. https://doi.org/10.1038/s41596-019-0285-9.
  Gryder, B. E. et al. Nat. Genet. 51, 1714–1722 (2019): https://doi.org/10.1038/s41588-019-0534-4 Stanton, B. et al. Preprint at https://protocolexchange.researchsquare.com/article/nprot-7121/v1 (2018): https://doi.org/10.1038/protex.2018.130
</details>

- [Pore-C tools](https://github.com/nanoporetech/pore-c/) - Pore-C chromatin conformation capture using Oxford Nanopore Technologies, direct sequencing of multi-way chromatin contacts without amplification (concatemers, HOLR - high-order and long-range contacts). >18 times higher efficiency as compared with SPRITE, enrichment in concatemers. Tested on the Gm12878 cell line. Hi-C matrix from Pore-C well resembles published data. Concatemers show significantly lower distance decay. Concatemers better resolve complex cancer rearrangements, well-suited for de novo genome scaffolding. [Snakemake pipeline](https://github.com/nanoporetech/Pore-C-Snakemake), [detection of multi-way interactions](https://github.com/mskilab/GxG). <details>
    <summary>Paper</summary>
    Ulahannan, Netha, Matthew Pendleton, Aditya Deshpande, Stefan Schwenk, Julie M. Behr, Xiaoguang Dai, Carly Tyer, et al. "Nanopore Sequencing of DNA Concatemers Reveals Higher-Order Features of Chromatin Structure"  https://doi.org/10.1101/833590 Preprint. Genomics, November 7, 2019. 
</details>

- [TiledC](https://github.com/oudelaar/TiledC) - Tiled-C low-input 3C technology, requiring >20,000 cells. Applied for in vivo mouse erythroid differentiation, alpha-globin genes. TADs are pre-existing, regulatory interactions gradually form during differentiation. Integration with scRNA-seq data (CITE-seq technology) and ATAC-seq data. Analyzed using [CCseqBasic pipeline](https://github.com/Hughes-Genome-Group/CCseqBasicF/). [Data (Tiled-C, CITE-seq, ATAC-seq)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137477). <details>
    <summary>Paper</summary>
    Oudelaar, A. Marieke, Robert A. Beagrie, Matthew Gosden, Sara de Ornellas, Emily Georgiades, Jon Kerry, Daniel Hidalgo, et al. "Dynamics of the 4D Genome during Lineage Specification, Differentiation and Maturation in Vivo"  https://doi.org/10.1101/763763 Preprint. Genomics, September 10, 2019. 
</details>

- [Red-C (RNA Ends on DNA Capture) technology](https://github.com/agalitsyna/RedClib) - proximity ligation-based technology to understand RNA-DNA interactome (3' and 5' ends of the RNA molecule associated with a given DNA site). Technology description(Methods, Results, Figure 1). Assesses all coding and noncoding RNA classes (piRNAs, tRNAs, rRNAs, snRNAs, scRNAs, srpRNAs, vlinc, antisense RNAs). Applied to K562 cells. The RNA-DNA matrices show a wide distribution of RNA contacts along extended genomic regions. The highest number of contacts was observed for mRNA, linc and vlinc RNAs. Enhancer RNAs (eRNAs) produced from strong enhancers interact with other strong enhancers on the same chromosome.  MIR3648 and MIR3687 interact with repressed chromatin genome-wide. Investigation of mRNA interaction preference across gene body. Discovery of novel chromatin-interacting RNAs (X-RNAs). [Processed Red-C and RNA-seq data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136141), [GitHub, Nextflow processing pipeline](https://github.com/agalitsyna/RedClib). <details>
    <summary>Paper</summary>
    Gavrilov, Alexey A, Anastasiya A Zharikova, Aleksandra A Galitsyna, Artem V Luzhin, Natalia M Rubanova, Arkadiy K Golov, Nadezhda V Petrova, et al. "Studying RNA–DNA Interactome by Red-C Identiﬁes Noncoding RNAs Associated with Various Chromatin Types and Reveals Transcription Dynamics"  https://doi.org/10.1093/nar/gkaa457 July 6, 2020
</details>

- [ChIA-Drop technology](https://doi.org/10.1038/s41586-019-0949-1) - multiplex chromatin interaction analysis via droplet-based and barcode-linked sequencing. 10X genomics technology, instead of cells, each gel-bead-in-emulsion contains crosslinked (but not ligated) chromatin and unique barcode. Singleton reads removed, distance constraints further refine multi-way interactions. FISH agrees with ChIA-Drop. Only approx. 20% of promoters are regulated by multiple enhancers, in Drosophila S2 cells. Data. Supplementary - notes about data processing, overview of [ChIA-Dropbox](https://github.com/TheJacksonLaboratory/ChIA-DropBox) software for visualization, Python implementation. [Drosophila S2 cell line ChIA-Drop data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109355), [Hi-C](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99104). <details>
    <summary>Paper</summary>
    Zheng, Meizhen, Simon Zhongyuan Tian, Daniel Capurso, Minji Kim, Rahul Maurya, Byoungkoo Lee, Emaly Piecuch, et al. "Multiplex Chromatin Interactions with Single-Molecule Precision"  https://doi.org/10.1038/s41586-019-0949-1 Nature (February 2019)
</details>

- [ChIA-Dropbox](https://github.com/TheJacksonLaboratory/ChIA-DropBox) - analysis and visualization of ChIA-Drop data (10X Genomics-based GEMs containing cross-linked chromatin). Alignment, demultiplexing, reconstruction of intra-chromosomal chromatin droplets from barcodes, drops inter-chromosomal, generates visualizations compatible with Juicebox, BASIC browser, its own ChIA-View (R/Shiny app). QC plots, multi-way interactions visualization, cluster/fragment/promoter views. Global and RNAPII-enriched [GEO GSE109355](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109355) - ChIA-Drop data for Drosophila. <details>
    <summary>Paper</summary>
    Zheng, Meizhen. "ChIA-DropBox: A Novel Analysis and Visualization Pipeline for Multiplex Chromatin Interactions"  https://doi.org/10.1101/613034 Preprint, April 18, 2019
</details>

- [dip-c - Tools to analyze Dip-C (or other 3C/Hi-C) data](https://github.com/tanlongzhi/dip-c) -  technology for single-cell Hi-C, multiplex end-tagging amplification (META). Detects approx. 5 times more contacts. Possible to make haplotype-separated Hi-C maps, detect CNVs, resolve X-chromosome inactivation. Approx. 10kb-resolution Hi-C matrices, 3D genome reconstruction at 20kb resolution. PCA on chromatin compartments separates cell types. Comparison with Nagano data, bulk Gm12878 Hi-C. [hickit](https://github.com/lh3/hickit). Processed data: [Gm12878, 17 cells, PBMC, 18 cells](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117876). <details>
    <summary>Paper</summary>
    Tan, Longzhi, Dong Xing, Chi-Han Chang, Heng Li, and X. Sunney Xie. "Three-Dimensional Genome Structures of Single Diploid Human Cells"  https://doi.org/10.1126/science.aat5641 Science. August 31, 2018
</details>

- [Methyl-HiC technology, in situ Hi-C and WGBS](https://bitbucket.org/dnaase/bisulfitehic/src/master/) - Comparable Hi-C matrices, TADs. 20% fewer CpGs overall, more CpGs in open chromatin. Proximal CpGs correlate irrespectively of loop anchors, weaker for inter-chromosomal interactions. Application to single-cell, mouse ESCs under different conditions. Relevant clustering, cluster-specific genes. Methods for wet-lab and computational processing. [Bulk (replicates) and single-cell Methyl-HiC data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE119171). Bhmem pipeline to map bisulfite-converted reads, Juicer pipeline for processing, VC normalization, HiCRep at 1Mb matrix similarity. <details>
    <summary>Paper</summary>
    Li, Guoqiang, Yaping Liu, Yanxiao Zhang, Naoki Kubo, Miao Yu, Rongxin Fang, Manolis Kellis, and Bing Ren. "Joint Profiling of DNA Methylation and Chromatin Architecture in Single Cells"  https://doi.org/10.1038/s41592-019-0502-z Nature Methods, August 5, 2019. 
</details>

- [Review of chromatin conformation capture technologies](https://doi.org/10.1038/s41576-019-0195-2) - from image-based methods (FISH), through proximity ligation (3/4/5C, Hi-C, TCC, ChIA-PET, scHi-C), to ligation-free methods (GAM, SPRITE, ChIA-Drop). Details of each technology (Table 1, Figures), comparison of them (Table 2). <details>
    <summary>Paper</summary>
    Kempfer, Rieke, and Ana Pombo. "Methods for Mapping 3D Chromosome Architecture"  https://doi.org/10.1038/s41576-019-0195-2 Nature Reviews Genetics, December 17, 2019. 
</details>

- [Optimization steps for Hi-C wet-lab protocol](https://doi.org/10.1101/287201) - Pitfalls and their effect on the downstream quality. Recommendations for each step. 
    Golloshi, Rosela, Jacob Sanders, and Rachel Patton McCord. "Iteratively Improving Hi-C Experiments One Step at a Time"  https://doi.org/10.1101/287201 Preprint. Genomics, March 22, 2018. 

- [DLO Hi-C technology (digestion-ligation-only Hi-C)](https://github.com/GangCaoLab/DLO-HiC-Tools) - Uses two rounds of digestion and ligation, without biotin and pull-down. Allows for early evaluation of Hi-C quality. Cost-effective, high signal-to-noise-ratio. Tested on THP-1 (human monocytes) and K562 cells. [Data processed with ChIA-PET Tool, normalized with ICE](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89663). <details>
    <summary>Paper</summary>
    Lin, Da, Ping Hong, Siheng Zhang, Weize Xu, Muhammad Jamal, Keji Yan, Yingying Lei, et al. "Digestion-Ligation-Only Hi-C Is an Efficient and Cost-Effective Method for Chromosome Conformation Capture"  https://doi.org/10.1038/s41588-018-0111-2 Nature Genetics 50, no. 5 (May 2018)
</details>

- [HIPMap](https://doi.org/10.1101/sqb.2015.80.027417)  - high-throughput imaging and analysis pipeline to map the location of gene loci within the 3D space. [FISH](#fish) in a 384-well plate format, automated imaging. <details>
    <summary>Paper</summary>
    Shachar, Sigal, Gianluca Pegoraro, and Tom Misteli. "HIPMap: A High-Throughput Imaging Method for Mapping Spatial Gene Positions"  https://doi.org/10.1101/sqb.2015.80.027417 Cold Spring Harbor Symposia on Quantitative Biology 80 (2015)
</details>

- [HiC method description](https://doi.org/10.1126/science.1181369) - 1Mb, Gm06990. Small chromosomes, but 18, interact. Compartment A associated with open chromatin. 1Mb, 100kb resolution. <details>
    <summary>Paper</summary>
    Lieberman-Aiden, Erez, Nynke L. van Berkum, Louise Williams, Maxim Imakaev, Tobias Ragoczy, Agnes Telling, Ido Amit, et al. "Comprehensive Mapping of Long-Range Interactions Reveals Folding Principles of the Human Genome"  https://doi.org/10.1126/science.1181369 Science (New York, N.Y.) 326, no. 5950 (October 9, 2009)
</details>

- [3C technology](https://doi.org/10.1126/science.1067799) - the matrix of interaction frequencies, application to reveal spatial information, applied to yeast (S cerevisiae) genome. Interphase and metaphase chromosomes show different patterns of interactions. Distance-dependent decay of interaction frequencies. Basic observations on chromosome size, inter-chromosomal interactions. <details>
    <summary>Paper</summary>
    Dekker, Job, Karsten Rippe, Martijn Dekker, and Nancy Kleckner. "Capturing Chromosome Conformation"  https://doi.org/10.1126/science.1067799 Science (New York, N.Y.) 295, no. 5558 (February 15, 2002)
</details>

- [Hi-C library prep with Arima kit](https://doi.org/10.1016/j.celrep.2021.109722) - The Arima kit uses two restriction enzymes: ^GATC (DpnII), G^ANTC (N can be either of the 4 genomic bases) (HinfI), which after ligation of DNA ends generates 4 possible ligation junctions in the chimeric reads: GATC-GATC, GANT-GATC, GANT-ANTC, GATC-ANTC. Source: 'Hi-C library preparation' Methods section from [Du et al., “DNA Methylation Is Required to Maintain Both DNA Replication Timing Precision and 3D Genome Organization Integrity.”"  https://doi.org/10.1016/j.celrep.2021.109722)

#### Micro-C

- Micro-Capture-C (MCC) technology, wet-lab and computational protocols. Introduction to the 3C technologies, advantages and limitations of each. MCC provides targeted data at base pair resolution, can be used to study promoters, enhancers, boundary elements, and help to interpret the effect of disease-associated variants in the noncoding genome. Requires 3M cells, 10^5-10^6 read depth. Oligo design using [Capsequm](https://github.com/jbkerry/capsequm), the full pipeline - [Micro Capture-C analysis tools](https://process.innovation.ox.ac.uk/software/p/16529a/micro-capture-c-academic/1), [GitHub](https://github.com/jojdavies/Micro-Capture-C). <details>
    <summary>Paper</summary>
    Hamley, Joseph C., Hangpeng Li, Nicholas Denny, Damien Downes, and James O. J. Davies. “Determining Chromatin Architecture with Micro Capture-C.” Nature Protocols, March 29, 2023. https://doi.org/10.1038/s41596-023-00817-8.
</details>

- Micro-Capture-C (Tiled-MCC, MNase digestion on viewpoints) technology to generate 3C interaction profiles at base-pair resolution. Mouse Embryonic Stem Cells. Local cis-regulatory loops depend on CTCF but not cohesin. Cohesin (RAD21) depletion slightly affects long-range interactions, CTCF depletion leads to rewiring of interactions. Loop extrusion is not essential for enhancer-promoter itteractions, but for their robustness and precision. [GitHub](https://github.com/jojdavies/Micro-Capture-C) and [the Oxford University Innovation Software Store](https://process.innovation.ox.ac.uk/software/p/16529a/micro-capture-c-academic/1) - scripts for analysing Micro Capture-C (MCC) data. [GSE181694](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181694) - Tiled MCC data, replicates of WT, untreated, RAD21-depleted and CTCF-depleted conditions. Plus, publicly available RNA-seq, CTCF/RAD21 ChIP-seq in [Data availability](https://www.nature.com/articles/s41467-022-29696-5#data-availability) section. [Tweet](https://twitter.com/MariekeOudelaar/status/1516342440227131392?s=20&t=I6UCnElJPmHGGrHsi3UWtw). <details>
    <summary>Paper</summary>
    Aljahani, Abrar, Peng Hua, Magdalena A. Karpinska, Kimberly Quililan, James OJ Davies, and A. Marieke Oudelaar. "Analysis of sub-kilobase chromatin topology reveals nano-scale regulatory interactions with variable dependence on cohesin and CTCF." Nature communications 13, no. 1 (2022): 1-13. https://doi.org/10.1038/s41467-022-29696-5 
</details>

- Stability of enhancer-promoter (E-P) interactions after 3-hours depletion of CTCF, cohesin, WAPL, and YY1 (auxin-inducible degron). Micro-C, nascent transcript profiling (mNET-seq), live-cell single-molecule imaging (single-particle tracking technique spaSPT. Mouse embryonic stem cells. E-P interaction strength correlates with gene expression, cohesin loop strength does not. E-P and P-P interactions can cross TAD boundaries. TADs are much more effective at insulating cohesin loops than E-P and P-P loops. CTCF depletion had minimal impact on genomewide interactions; RAD21 depletion reduced contact frequency in the range of 10 – 200 kb but increased interactions at 300 kb – 5 Mb; and WAPL depletion increased contacts at 70 – 700 kb but reduced contacts at 1 – 5 Mb. YY1 depletion has minimal effect on gene expression and E-R and P-P interactions. [GEO GSE130257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130275) - Previously published Micro-C replicates of mESCs WT and RNA PolII inhibition, .mcool and .hic files. [Zenodo 10.5281/zenodo.5035837](https://zenodo.org/record/5035837) - spaSPT data for endogenously tagging HaloTag-YY1 in mESCs. [GEO GSE178982](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178982) - ChIP-seq, RNA-seq, mNET-seq, Micro-C replicates of CTCF, RAD21, WAPL, YY1 depletion, pooled .mcool files. <details>
    <summary>Paper</summary>
    Hsieh, Tsung-Han S., Claudia Cattoglio, Elena Slobodyanyuk, Anders S. Hansen, Xavier Darzacq, and Robert Tjian. “Enhancer-Promoter Interactions and Transcription Are Maintained upon Acute Loss of CTCF, Cohesin, WAPL, and YY1.” Preprint. Molecular Biology, July 14, 2021. https://doi.org/10.1101/2021.07.14.452365.
</details>

- Micro-Capture-C (MCC) technology. Advanced and optimized 3C technology, using micrococcal nuclease (MNase), intact cells permeabilized with digitonin, ultra-deep sequencing, plus sequencing the ligation junction. 330 viewpoints (promoters, enhancers, CTCF binding sites) in primary mouse erythroid and empryonic stem cells. Many biological observations, cohesin is required for CTCF contacts.Enhancer-promoter contacts are mediated by at least partly different mechanisms than cohesin-mediated interactions. [Computational pipeline](https://process.innovation.ox.ac.uk/software/p/16529a/micro-capture-c-academic/1). [GEO GSE153256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153256) - mouse spleen erythroid cells MCC data. [GEO GSE144336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144336) - mouse primary erythroid cells, embryonic stem cells, HUDEP-2 cells. MCC, CTCF/Rad21/NIPBL ChIP-seq. <details>
    <summary>Paper</summary>
    Hua, Peng, Mohsin Badat, Lars L. P. Hanssen, Lance D. Hentges, Nicholas Crump, Damien J. Downes, Danuta M. Jeziorska, et al. "Defining Genome Architecture at Base-Pair Resolution"  https://doi.org/10.1038/s41586-021-03639-4 Nature, (July 1, 2021)
</details>

- Ultra-deep Micro-C maps of human embryonic stem cells and fibroblasts. Compared with in situ Hi-C (DpnII 4bp cutter), Micro-C allows for the detection of \~20,000 additional loops, improved signal-to-noise ratio. Similar distance-dependent decay, recovery of A/B compartments, better recovery of close-range interactions. High-resolution interaction boundaries are not created equal - most are CTCF+ and YY1+, some are CTCF- and YY1+, CTCF- and YY1-, and completely negative boundaries. Multiple, weak pause sites of SMC complexes. Distiller-nf processing, ICE normalization, other Mirnylab tools. [Data](https://data.4dnucleome.org/). <details>
    <summary>Paper</summary>
    Krietenstein, Nils, Sameer Abraham, Sergey V. Venev, Nezar Abdennur, Johan Gibcus, Tsung-Han S. Hsieh, Krishna Mohan Parsi, et al. "Ultrastructural Details of Mammalian Chromosome Architecture"  https://doi.org/10.1016/j.molcel.2020.03.003 Molecular Cell 78, no. 3 (May 2020)
</details>

- Micro-C (MNase digestion Hi-C) data analysis. Nucleosome-level (\~100-200bp) resolution Hi-C, captures all structures in regular Hi-C data. Stripes/flames structures correspond to enhancer-promoter interactions, colocalize with PolII, CTCF.  TADs are further split into "micro TADs", insulation score. Active boundaries colocalize with CpG islands, promoters, tRNA genes. Inactive boundaries are in repeats. Micro TADs subdivided into five subgroups. Two-start zig-zag model tetra-nucleosome stacks. Mouse stem cell, 38 biological replicates. [Twitter](https://twitter.com/Stanley_TH/status/1242843692081111040?s=20). <details>
    <summary>Paper</summary>
    Hsieh, Tsung-Han S., Claudia Cattoglio, Elena Slobodyanyuk, Anders S. Hansen, Oliver J. Rando, Robert Tjian, and Xavier Darzacq. "Resolving the 3D Landscape of Transcription-Linked Mammalian Chromatin Folding"  https://doi.org/10.1016/j.molcel.2020.03.002 Molecular Cell, March 2020
</details>

- Micro-C (MNase digestion Hi-C) technology and basic analysis. Human embryonic stem cells H1-ESC and differentiated human foreskin fibroblasts (HFFc6). Captures standard Hi-C features, with many additional interaction peaks ("dots"). Enrichment of classical marks of TAD boundaries (Fig 3C) - RAD21, TAF1, PHF8, CTCF, TBP, POL2RA, YY1, and more. <details>
    <summary>Paper</summary>
    Krietenstein, Nils, Sameer Abraham, Sergey Venev, Nezar Abdennur, Johan Gibcus, Tsung-Han Hsieh, Krishna Mohan Parsi, et al. "Ultrastructural Details of Mammalian Chromosome Architecture"  https://doi.org/10.1101/639922 Preprint. Genomics, May 17, 2019. 
</details>

- Micro-C technology - mononucleosome resolution mapping in yeast. Micrococcal nuclease to fragment chromatin. Yeast does not have TADs, but Micro-C revealed self-associating domains (chromatin interaction domains, CIDs) driven by the number of genes. Enrichment of histone modifications. [Data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68016). <details>
    <summary>Paper</summary>
    Hsieh, Tsung-Han S., Assaf Weiner, Bryan Lajoie, Job Dekker, Nir Friedman, and Oliver J. Rando. "Mapping Nucleosome Resolution Chromosome Folding in Yeast by Micro-C"  https://doi.org/10.1016/j.cell.2015.05.048 Cell 162, no. 1 (July 2015)
</details>

- [Micro-C](https://github.com/dovetail-genomics/Micro-C) scripts and [documentation](https://micro-c.readthedocs.io/) for Dovetail™ Micro-C Kit. Uses the Micrococcal nuclease (MNase) enzyme instead of restriction enzymes for chromatin digestion, yielding 146 bp fragments distributed frequently across the genome. Detailed computational processing of Micro-C data, generation of valid pairs using bwa mem, pairtools, juicer_tools, cooler and other genomics tools, QC, contact matrix generation. [Introduction to Dovetail™ Micro-C Technology Video](https://vimeo.com/431456105), 5min, [Breaking the Hi-C Resolution Barrier with Micro-C Video](https://vimeo.com/470724936), 45 min. [Micro-C datasets](https://micro-c.readthedocs.io/en/latest/data_sets.html), [MicroC](https://github.com/EricSDavis/MicroC) Snakemake pipeline by Eric Davis.


#### Multi-way interactions

- **ImmunoGAM** - an extension of genome architecture mappint (GAM, thin nuclear cryosections. Works with small numbers of cells (approx. 400-1500 cells). Immunoselection of cryosections (incubation with antibodies). Applied to profile oligodendrocytes (OLGs), pyramidal glutamatergic neurons (PGNs), dopaminergic neurons (DNs) of mouse brain. Integration with pseudobulk RNA-seq and ATAC-seq. Insulation score to detect domains, extensive cell type-specific TAD, AB compartment rearrangement. "Melting" - contact loss at highly expressed long neuronal genes (MELTRON pipeline to calculate a melting score https://github.com/DominikSzabo1/Meltron). Methods: GAM data window calling and QC. [Extensive supplementary material](https://www.nature.com/articles/s41586-021-04081-2#Sec55). [GEO GSE148792](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148792) - GAM sequencing data from DN, PGN and OLG, non-normalized co-segregation matrices, normalized pair-wised chromatin contacts maps and raw GAM segregation tables. [Microscopy images](https://github.com/pombo-lab/WinickNg_Kukalev_Harabula_Nature_2021/tree/main/microscopy_images/). A public UCSC session with all data produced, as well as all published data utilized in this study is available at [the UCSC Genome Browser session](http://genome-euro.ucsc.edu/s/Kjmorris/Winick_Ng_2021_GAMbrainpublicsession). [Tweet](https://twitter.com/LabPombo/status/1461008979837005826?s=20). <details>
    <summary>Paper</summary>
    Winick-Ng, Warren, Alexander Kukalev, Izabela Harabula, Luna Zea-Redondo, Dominik Szabó, Mandy Meijer, Leonid Serebreni, et al. "Cell-Type Specialization Is Encoded by Specific Chromatin Topologies"  https://doi.org/10.1038/s41586-021-04081-2 Nature, November 17, 2021. 
</details>

- [MATCHA](https://github.com/ma-compbio/MATCHA) - a computational framework for the analysis of multi-way chromatin interactions (hyperedges). Graph embedding fed into Hyper-SAGNN (self-attention based graph neural network for hypergraphs). Applied to SPRITE and ChIA-Drop data, Gm12878 at 1Mb and 100kb resolution. <details>
    <summary>Paper</summary>
    Zhang, Ruochi, and Jian Ma. "MATCHA: Probing Multi-Way Chromatin Interaction with Hypergraph Representation Learning"  https://doi.org/10.1016/j.cels.2020.04.004 Cell Systems, (May 2020)
</details>

- **Multi-contact 3C** (MC-3C, based on conventional 3C, C-walk, and multi-contact 4C approaches) technology reveals distinct chromosome territories with very little mixing, never entanglement, same with chromosomal compartment domains (A-A, B-B interactions predominant, A-B - minimal to none). Analysis of C-walks - connected paths of pairwise interactions. Compared with C-walks generated from Hi-C data. Permutation analysis of the significance of insulated, mixed, and intermediate domains. PacBio sequencing, processing with SMRT Analysis package and a [custom pipeline](https://github.com/dekkerlab/MC-3C_scripts). <details>
    <summary>Paper</summary>
    Tavares-Cadete, Filipe, Davood Norouzi, Bastiaan Dekker, Yu Liu, and Job Dekker. "Multi-Contact 3C Data Reveal That the Human Genome Is Largely Unentangled"  https://doi.org/10.1101/2020.03.03.975425 Preprint. Genomics, March 4, 2020. 
</details>

- [MIA-Sig](https://github.com/TheJacksonLaboratory/mia-sig) (multiplex interactions analysis by signal processing algorithms) - de-noise the multi-interaction data, assess the statistical significance of chromatin complexes, identify topological domains and multi-way inter-TAD contacts. Tailored for ChIA-Drop, also works with SPRITE and GAM data. Uses the assumption that the genomic distance of a true multiway interaction should be more evenly spaced to remove noise. <details>
    <summary>Paper</summary>
    Kim, Minji, Meizhen Zheng, Simon Zhongyuan Tian, Byoungkoo Lee, Jeffrey H. Chuang, and Yijun Ruan. "MIA-Sig: Multiplex Chromatin Interaction Analysis by Signalprocessing and Statistical Algorithms"  https://doi.org/10.1186/s13059-019-1868-z Genome Biology, (December 2019)
</details>

- **MC-4C** multi-way contacts technology and computational protocols. ~2 weeks, ~$600/sample, best for <120kb regions. [Computational protocol](https://github.com/deLaatLab/mc4c_py), test data included. <details>
    <summary>Paper</summary>
    Vermeulen, Carlo, Amin Allahyar, Britta A. M. Bouwman, Peter H. L. Krijger, Marjon J. A. M. Verstegen, Geert Geeven, Christian Valdes-Quezada, et al. "Multi-Contact 4C: Long-Molecule Sequencing of Complex Proximity Ligation Products to Uncover Local Cooperative and Competitive Chromatin Topologies"  https://doi.org/10.1038/s41596-019-0242-7 Nature Protocols 15, no. 2 (February 2020)
</details>

- **MC-4C** - Multi-way interactions technology, uses Nanopore MinION (or, PacBio) sequencing. Cross-linking, cutting with four-cutter and six-cutter enzymes, circularization, cutting with Cas9 gRNA designed to the viewpoint region, selective amplification of concatemers with primers specific to the viewpoint. Rigorous filtering strategy, interactions are allowing to distinguish reads coming from individual alleles. Compared with genome-wide multi-contact technologies C-walks, SPRITE, GAM, Tri-C. Applied to mouse beta-globin (fetal liver where hemoglobin genes are expressed, and brain, where they are silent) and protocadherin-alpha (same tissues, vice versa ) loci. Super enhancers can form hubs, target multiple genes. WAPL deletion in HAP1 (leukemia) cells stimulates the collision of CTCF-anchored domain loops to form rosette-like structures. [MC-4C processing pipeline](https://github.com/UMCUGenetics/pymc4c/), [Visualization of the analyzed data](http://www.multicontactchromatin.nl/), [raw data](https://www.ebi.ac.uk/ena/data/view/PRJEB23327), [processed data matrices](https://data.mendeley.com/datasets/wbk8hk87r2/1). <details>
    <summary>Paper</summary>
    Allahyar, Amin, Carlo Vermeulen, Britta A. M. Bouwman, Peter H. L. Krijger, Marjon J. A. M. Verstegen, Geert Geeven, Melissa van Kranenburg, et al. "Enhancer Hubs and Loop Collisions Identified from Single-Allele Topologies"  https://doi.org/10.1038/s41588-018-0161-5 Nature Genetics 50, no. 8 (August 2018)
</details>

- **Multiplex-GAM** - Genome Architecture Mapping technology. Multiple nuclear profiles (2-3 NPs, slices of nuclei) are sequenced together. Updates of the SLICE method for GAM data analysis. Details on the dependencies between the number of GAM samples, probability to detect interactions, the number of NPs per sample, technical parameters. TAD boundaries match those from Hi-C. Method-specific contacts exist, but GAM contacts are enriched in CTCF-enhancer contacts and other active transcription annotations, in contrast to heterochromatin-prone Hi-C contacts. No data yet. <details>
    <summary>Paper</summary>
    Beagrie, Robert A, Christoph J Thieme, Carlo Annunziatella, Catherine Baugher, Yingnan Zhang, Markus Schueler, Dorothee CA Kramer, et al. "Multiplex-GAM: Genome-Wide Identification of Chromatin Contacts Yields Insights Not Captured by Hi-C"  https://doi.org/10.1101/2020.07.31.230284 Preprint. Molecular Biology, July 31, 2020. 
</details>

- **GAM** (genome architecture mapping) - restriction- and ligation free chromatin conformation capture technology. Isolates and sequences DNA content of many ultra-thin (~0.22um) cryo-fixes nuclear slices, whole-genome amplification.~30kb matrix is built on frequencies of co-occurrences of regions in multiple slices. Applied to mouse ESCs. Normalization using linkage disequilibrium (better than ICE). GAM and Hi-C matrices are highly correlated, A/B compartments, TADs significantly overlap. Enhancers and active genes significantly interact. Multi-way interactions, triplets, super-enhancers are enriched in three-way interactions, confirmed by FISH. SLICE (Statistical Inference of co-segregation) mathematical model to assign significance of interactions ([Supplementary Note 1](https://static-content.springer.com/esm/art%3A10.1038%2Fnature21411/MediaObjects/41586_2017_BFnature21411_MOESM22_ESM.pdf)). Negative binomial modeling of significant interactions, log-normal noise. [Supplementary Table 2](https://static-content.springer.com/esm/art%3A10.1038%2Fnature21411/MediaObjects/41586_2017_BFnature21411_MOESM24_ESM.xlsx) - genomic coordinates of triplet (three-way interacting) TADs. Raw data, processed matrices, [Python scripts](https://gam.tools/papers/nature-2017/), [GEO GSE64881](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64881) - mouse ES cell GAM matrices at 1Mb and 30kb resolution. <details>
    <summary>Paper</summary>
    Beagrie, Robert A., Antonio Scialdone, Markus Schueler, Dorothee C. A. Kraemer, Mita Chotalia, Sheila Q. Xie, Mariano Barbieri, et al. "Complex Multi-Enhancer Contacts Captured by Genome Architecture Mapping"  https://doi.org/10.1038/nature21411 Nature 543, no. 7646 (23 2017)
</details>

- **SPRITE** (Split-Pool Recognition of Interactions by Tag Extension) - Multi-way interactions technology, barcode system based on repeated pooling ans splitting of crosslinked DNAhairballs to uncover clustered DNA fragments (SPRITE clusters). Does not involve proximity ligation. Transcriptionally active interaction hubs around duclear speckles, inactive - around nucleolus. These nuclear bodies constrain the 3D packaging. Applied to mouse embryonic stem cells (mESCs) and GM12878, data correlates well with Hi-C. [Data: 25kb-1Mb inter- and intra-chromosomal interactions](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114242), [Processing Pipeline](https://github.com/GuttmanLab/sprite-pipeline/wiki). <details>
    <summary>Paper</summary>
    Quinodoz, Sofia A., Noah Ollikainen, Barbara Tabak, Ali Palla, Jan Marten Schmidt, Elizabeth Detmar, Mason M. Lai, et al. "Higher-Order Inter-Chromosomal Hubs Shape 3D Genome Organization in the Nucleus"  https://doi.org/10.1016/j.cell.2018.05.024 Cell 174, no. 3 (July 2018)
</details>

- **C-walks**, multi-way technology, genome-wide. TADs organize chromosomal territories. Active and inactive TAD properties. Methods: Good mathematical description of insulation score calculations. Filter TADs smaller than 250kb. Inter-chromosomal contacts are rare, \~7-10%. Concatemers (more than two contacts) are unlikely. <details>
    <summary>Paper</summary>
    Olivares-Chauvet, Pedro, Zohar Mukamel, Aviezer Lifshitz, Omer Schwartzman, Noa Oded Elkayam, Yaniv Lubling, Gintaras Deikus, Robert P. Sebra, and Amos Tanay. "Capturing Pairwise and Multi-Way Chromosomal Conformations Using Chromosomal Walks"  https://doi.org/10.1038/nature20158 Nature 540, no. 7632 (November 30, 2016)
</details>

- **TM3C** - Tethered multiple 3C technology to probe multi-point contacts. NHEK, KBM7 cells, detected the Philadelphia chromosome, investigated triple contacts in the IGF2-H19 locus at 40kb, detected typical genomic structures (chromosomal compartments, distance-decay, TADs), reconstructed 3D genome at 1Mb resolution. A two-phase mapping strategy that separately maps chimeric subsequences within a single read (Methods). Multiple 4-cutter restriction enzymes. <details>
    <summary>Paper</summary>
    Ay, Ferhat, Thanh H Vu, Michael J Zeitz, Nelle Varoquaux, Jan E Carette, Jean-Philippe Vert, Andrew R Hoffman, and William S Noble. "Identifying Multi-Locus Chromatin Contacts in Human Cells Using Tethered Multiple 3C"  https://doi.org/10.1186/s12864-015-1236-7 BMC Genomics 16, no. 1 (2015)
</details>

- **INGRID** - chromatin conformation capture technology using in-gel replication of interacting DNA segments. Detects multi-way interactions. Protocol, demonstration of beta-globin gene locus. <details>
    <summary>Paper</summary>
    Gavrilov, Alexey A., Helena V. Chetverina, Elina S. Chermnykh, Sergey V. Razin, and Alexander B. Chetverin. "Quantitative Analysis of Genomic Element Interactions by Molecular Colony Technique"  https://doi.org/10.1093/nar/gkt1322 Nucleic Acids Research 42, no. 5 (March 1, 2014)
</details>

#### Imaging

- 3D chromatin imaging technology review. [Table 1](https://www.cell.com/action/showFullTableHTML?isHtml=true&tableId=tbl1&pii=S2211-1247%2823%2900372-8)- 3D technologies and their description. Details in the text. Live-cell imaging, spatial and temporal chromatin profiling. <details>
  <summary>Paper</summary>
  Cosma, Maria Pia, and Maria Victoria Neguembor. “The Magic of Unraveling Genome Architecture and Function.” Cell Reports, April 2023, 112361. https://doi.org/10.1016/j.celrep.2023.112361.
</details>

- [OligoFISSEQ](https://doi.org/10.1038/s41592-020-0890-0) - fluorescence in situ sequencing (FISSEQ) of barcoded oligopaint probes to enable the rapid visualization of many targeted genomic regions. Three strategies for interrogation of targets: sequencing by ligation (SBL, O-LIT), sequencing by synthesis (SBS, O-SIT), sequencing by hybridization (SBH, O-HIT). Applied to human diploit fibroblast cells, 36 regions across 6 chromosomes in just 4-8 rounds of sequencing. Fine tracing of 46 regions on X chromosome. Combined with OligoSTORM, allows for single-molecule super-resolution imaging. Whole genome-ready technology. [Figure-specific source data](https://www.nature.com/articles/s41592-020-0890-0#Sec51), [Analysis scripts](https://github.com/3DGenomes/OligoFISSEQ/). <details>
    <summary>Paper</summary>
    Nguyen, Huy Q., Shyamtanu Chattoraj, David Castillo, Son C. Nguyen, Guy Nir, Antonios Lioutas, Elliot A. Hershberg, et al. "3D Mapping and Accelerated Super-Resolution Imaging of the Human Genome Using in Situ Sequencing"  https://doi.org/10.1038/s41592-020-0890-0 Nature Methods, (August 2020)
</details>

- 3D genome visualization techniques overview. Electron microscopy, STORM and its modifications (Oligopaint and STORM = OligoSTORM). Nucleosome clutches and clutch domains, correlate with histone acetylation. A/B compartments and TADs generally correlate with imaging, but heterogeneity at single-cell level is observed. Live-cell imaging techniques (LacO/LacI, ParB/parS, CRISPR-based). Newer high-resolution techniques (3D MINFLUX). <details>
  <summary>Paper</summary>
  Lakadamyali, Melike, and Maria Pia Cosma. “Visualizing the Genome in High Resolution Challenges Our Textbook Understanding.” Nature Methods, March 2, 2020. https://doi.org/10.1038/s41592-020-0758-3.
</details>

- [MERFISH](https://github.com/BogdanBintu/ChromatinImaging) - Super-resolution imaging technology, reconstruction 3D structure in single cells at 30kb resolution, 1.2Mb region of Chr21 in IMR90 cells. Distance maps obtained by microscopy show small distance for loci within, and larger between, TADs. TAD-like structures exist in single cells. 2.5Mb region of Chr21 in HCT116 cells, cohesin depletion does not abolish TADs, only alter their preferential positioning. Multi-point (triplet) interactions are prevalent. TAD boundaries are highly heterogeneous in single cells. , diffraction-limited and STORM (stochastic optical reconstruction microscopy) imaging. [GitHub](https://github.com/BogdanBintu/ChromatinImaging). <details>
    <summary>Paper</summary>
    Bintu, Bogdan, Leslie J. Mateo, Jun-Han Su, Nicholas A. Sinnott-Armstrong, Mirae Parker, Seon Kinrot, Kei Yamaya, Alistair N. Boettiger, and Xiaowei Zhuang. "Super-Resolution Chromatin Tracing Reveals Domains and Cooperative Interactions in Single Cells"  https://doi.org/10.1126/science.aau1783 Science, (October 26, 2018)
</details>

### Normalization

- Lyu, Hongqiang, Erhu Liu, and Zhifang Wu. "[Comparison of Normalization Methods for Hi-C Data](https://europepmc.org/article/med/31588782) BioTechniques 68, no. 2 (2020) - a comprehensive analysis of six Hi-C normalization methods for their ability to remove systematic biases. The introduction provides a good classification and overview of different normalization methods, including the latest methods for cross-sample normalization, such as [multiHiCcompare](#multihiccompare). Human and mouse Hi-C data were used, only cis interaction matrices are considered. A systematic protocol for benchmarking is presented. Several benchmarks were performed, including statistical quality, the influence of resolution, consistency of distance-dependent changes in interaction frequency, reproducibility of the 3D architecture. [multiHiCcompare](#multihiccompare) is reported as outperforming other methods on a range of performance metrics. 

- Imakaev, Maxim, Geoffrey Fudenberg, Rachel Patton McCord, Natalia Naumova, Anton Goloborodko, Bryan R. Lajoie, Job Dekker, and Leonid A. Mirny. "[Iterative Correction of Hi-C Data Reveals Hallmarks of Chromosome Organization](https://doi.org/10.1038/nmeth.2148)" Nature Methods 9, no. 10 (October 2012) - ICE - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. Assumption - all loci should have equal visibility. Deconvolution into eigenvectors/values.

- Cournac, Axel, Hervé Marie-Nelly, Martial Marbouty, Romain Koszul, and Julien Mozziconacci. “[Normalization of a Chromosomal Contact Map](https://doi.org/10.1186/1471-2164-13-436).” BMC Genomics 13 (2012) - Normalization (Sequential Component Normalization, SCN) of Hi-C matrices. Three types of biases (GC content, restriction fratment length, new bias from circularization of DNA molecules). Technology overview, types of ligation fragments (Figure 1). 7% inter- and 14% intrachromosomal long-range interactions recovery, >80% removed. Normalization enhances the detection of CTCF-bound interactions. SCN algorithm (Methods) and implementation [HiCcompare::SCN](https://rdrr.io/bioc/HiCcompare/man/SCN.html) by [John Stansfield](https://github.com/jstansfield0).

- Yaffe, Eitan, and Amos Tanay. "[Probabilistic Modeling of Hi-C Contact Maps Eliminates Systematic Biases to Characterize Global Chromosomal Architecture](https://doi.org/10.1038/ng.947)" Nature Genetics 43, no. 11 (November 2011) - Sources of biases: 1) non-specific ligation (large distance between pairs); 2) length of each ligated fragments; 3) CG content and nucleotide composition; 4) Mappability. Normalization. Enrichment of long-range interactions in active promoters. General aggregation of active chromosomal domains. Chromosomal territories, high-activity and two low-activity genomic clusters

### Spectral clustering

- Y. X Rachel Wang, Purnamrita Sarkar, Oana Ursu, Anshul Kundaje and Peter J. Bickel, "Network modelling of topological domains using Hi-C data](https://arxiv.org/abs/1707.09587)"  - TAD analysis using graph-theoretical (network-based) methods. Treats TADs as a "community" within the network. Shows that naive spectral clustering is generally ineffective, leaving gaps in the data. 

- Liu, Sijia, Pin-Yu Chen, Alfred Hero, and Indika Rajapakse. "Dynamic Network Analysis of the 4D Nucleome"  https://doi.org/10.1101/268318 BioRxiv, January 1, 2018. - Temporal Hi-C data analysis using graph theory. Integrated with RNA-seq data. Network-based approaches such as von Neumann graph entropy, network centrality, and multilayer network theory are applied to reveal universal patterns of the dynamic genome. Toeplitz normalization. Graph Laplacian matrix. Detailed statistics.

- Norton, Heidi K, Harvey Huang, Daniel J Emerson, Jesi Kim, Shi Gu, Danielle S Bassett, and Jennifer E Phillips-Cremins. "Detecting Hierarchical 3-D Genome Domain Reconfiguration with Network Modularity"  https://doi.org/10.1101/089011 November 22, 2016. - Graph theory for TAD identification. Louvain-like local greedy algorithm to maximize network modularity. Vary resolution parameter, hierarchical TAD identification. Hierarchical spatial variance minimization method. ROC analysis to quantify performance. Adjusted RAND score to quantify the TAD overlap.

- Chen, Jie, Alfred O. Hero, and Indika Rajapakse. "Spectral Identification of Topological Domains"  https://doi.org/10.1093/bioinformatics/btw221 Bioinformatics (Oxford, England) 32, no. 14 (15 2016) - Spectral algorithm to define TADs. Laplacian graph segmentation using Fiedler vector iteratively. Toeplitz normalization to remove distance effect. Spectral TADs do not overlap with Dixon's, but better overlap with CTCF. Python implementation https://github.com/shappiron/TAD-Laplacian-Identification

## Courses

- [3DGenomes/3DAROC](https://github.com/3DGenomes/3DAROC) - list of notebooks for 4 days course on Hi-C data handling and 3D modeling of chromatin. Jupiter notebooks, using `pytadbit`. Instructors: Marc Marti-Renom, Francois Serra , David Castillo.

- [3d-genome-processing-tutorial](https://github.com/hms-dbmi/3d-genome-processing-tutorial) - A 3D genome data processing tutorial for ISMB/ECCB 2017. https://github.com/hms-dbmi/3d-genome-processing-tutorial

- [Workshop on measuring, analyzing, and visualizing the 3D genome with Hi-C data](https://github.com/hms-dbmi/hic-data-analysis-bootcamp). Presentations (PDFs, PPTX) and Jupyter notebooks. Cooler, HiCGlass, HiPlier. https://github.com/hms-dbmi/hic-data-analysis-bootcamp

- [The 3D Organization of Our Genome](https://youtu.be/Pl44JjA--2k) 3 min video recapitulates our current understanding of genome organization in the three-dimensional space of the cell nucleus. More videos in the description. By the [Cavalli lab videos](https://www.youtube.com/channel/UCFjv-tOfb_AaGLE-tKvLvTA), [Tweet](https://twitter.com/giacomo_cavalli/status/1453771881740476422?s=20)

## Labs

- [Best-Labs-of-3D-Genome](https://github.com/XiaoTaoWang/Best-Labs-of-3D-Genome) - Most active 3D genomics labs, collaborative networks. By [Xiaotao Wang](https://github.com/XiaoTaoWang)

<center><img src="https://github.com/XiaoTaoWang/Best-Labs-of-3D-Genome/blob/master/networks/2017-2021/2017-2021.coauthors.png?raw=true" alt="Best-Labs-of-3D-Genome" width="514" height="300"></center>

## Misc

- [PaintSHOP](https://paintshop.io/) - database and web tool for developing probes for oligonucleotide (oligo)-based FISH experiments. Uses OligoMiner, wrapped in Snakemake. RNA- and DNA-probe design. Newly generated collections for human, mouse, and model organisms. [GitHub1](https://github.com/beliveau-lab/PaintSHOP), [GitHub2](https://github.com/beliveau-lab/PaintSHOP_pipeline). <details>
  <summary>Paper</summary>
  Hershberg, Elliot A., Conor K. Camplisson, Jennie L. Close, Sahar Attar, Ryan Chern, Yuzhen Liu, Shreeram Akilesh, Philip R. Nicovich, and Brian J. Beliveau. “PaintSHOP Enables the Interactive Design of Transcriptome- and Genome-Scale Oligonucleotide FISH Experiments.” Nature Methods 18, no. 8 (August 2021): 937–44. https://doi.org/10.1038/s41592-021-01187-3.
</details>

- [HiC Poweraid](http://phanstiel-lab.med.unc.edu/poweraid/) - power analysis for loop detection from Hi-C data. Power calculated using the median counts for each genomic distance as depth values. Loop size is anti-correlated with fold change compression. Web app to assess power across sequencing depth and loop sizes. At least 6 billion valid contacts per condition, split between two replicates is required. <details>
  <summary>Paper</summary>
  Parker, Sarah M, Eric S Davis, and Douglas H Phanstiel. “Guiding the Design of Well-Powered Hi-C Experiments to Detect Differential Loops.” Preprint. Bioinformatics, March 16, 2023. https://doi.org/10.1101/2023.03.15.532762.
</details>

- [pairLiftOver](https://github.com/XiaoTaoWang/pairLiftOver) - Python package that converts the two-dimensional genomic coordinates of chromatin contact matrices (pairs) between genomic assemblies. Supports cool, [Juicer's](#juicer) hic, [HiC-Pro's](#hic-pro) allValidPairs, 4DN's pairs formats. Lifted matrices nearly equivalent realigned. By [Xiaotao Wang](https://github.com/XiaoTaoWang)

- [liftOverBedpe](https://github.com/dphansti/liftOverBedpe) - A liftOver wrapper to convert BEDPE files. Requires Python 2, e.g., `conda create -n liftOverBedpe python=2.7 -y`. By [Doug Phanstiel](https://github.com/dphansti)
