## Pipelines for Hi-C data processing

- `distiller-nf` - a modular Hi-C mapping pipeline for reproducible data analysis, nextflow version. Alignment, filtering, aggregating Hi-C matrices. https://github.com/mirnylab/distiller-nf

- `HiCExplorer` - set of programs to process, normalize, analyze and visualize Hi-C data, Python. https://hicexplorer.readthedocs.io/en/latest/, https://github.com/deeptools/HiCExplorer/

- `hiclib` - tools to qc, map, normalize, filter and analyze Hi-C data, https://bitbucket.org/mirnylab/hiclib

- `HiC-Pro` - An optimized and flexible pipeline for Hi-C data processing, https://github.com/nservant/HiC-Pro
    - Servant, Nicolas, Nelle Varoquaux, Bryan R. Lajoie, Eric Viara, Chong-Jian Chen, Jean-Philippe Vert, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiC-Pro: An Optimized and Flexible Pipeline for Hi-C Data Processing.” Genome Biology 16 (December 1, 2015): 259. https://doi.org/10.1186/s13059-015-0831-x. - HiC pipeline, references to other pipelines, comparison. From raw reads to normalized matrices. Normalization methods, fast and memory-efficient implementation of iterative correction normalization (ICE). Data format. Using genotyping information to phase contact maps.

- `HiCUP` pipeline, alignment only, output - BAM files. http://www.bioinformatics.babraham.ac.uk/projects/hicup/
    - Wingett, Steven, Philip Ewels, Mayra Furlan-Magaril, Takashi Nagano, Stefan Schoenfelder, Peter Fraser, and Simon Andrews. “HiCUP: Pipeline for Mapping and Processing Hi-C Data.” F1000Research 4 (2015): 1310. https://doi.org/10.12688/f1000research.7334.1. - Details about Hi-C sequencing artefacts. Used in conjunction with other pipelines.

- `HiTC` - High Throughput Chromosome Conformation Capture analysis, https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. https://bioconductor.org/packages/release/bioc/html/HiTC.html and https://www.bioconductor.org/packages/devel/bioc/vignettes/HiTC/inst/doc/HiTC.pdf 

- `my5C` tools - well-documented analysis and visualization of 5S data, http://my5c.umassmed.edu/

