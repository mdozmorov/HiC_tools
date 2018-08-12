## Pipelines for Hi-C data processing

- `distiller-nf` - Java modular Hi-C mapping pipeline for reproducible data analysis, nextflow pipeline. Alignment, filtering, aggregating Hi-C matrices. https://github.com/mirnylab/distiller-nf

- `Juicer` - Java full pipeline to convert raw reads into Hi-C maps, visualized in Juicebox. Call domains, loops, CTCF binding sites. `.hic` file format for storing multi-resolution Hi-C data. https://github.com/theaidenlab/juicebox/wiki/Download
    - Durand, Neva C., Muhammad S. Shamim, Ido Machol, Suhas S. P. Rao, Miriam H. Huntley, Eric S. Lander, and Erez Lieberman Aiden. “Juicer Provides a One-Click System for Analyzing Loop-Resolution Hi-C Experiments.” Cell Systems 3, no. 1 (July 2016): 95–98. https://doi.org/10.1016/j.cels.2016.07.002.
    - Rao, Suhas S. P., Miriam H. Huntley, Neva C. Durand, Elena K. Stamenova, Ivan D. Bochkov, James T. Robinson, Adrian L. Sanborn, et al. “A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping.” Cell 159, no. 7 (December 18, 2014): 1665–80. https://doi.org/10.1016/j.cell.2014.11.021. - Juicer analysis example. TADs defined by frequent interactions. Enriched in CTCF and cohesin members. Five domain types. A1 and A2 enriched in genes. Chr 19 contains 6th pattern B6. Enrichment in different histone modification marks. TADs are preserved across cell types. Yet, differences between Gm12878 and IMR90 were detected. Boundaries detection by scanning image. Refs to the original paper.


- `HiCExplorer` - set of Python scripts to process, normalize, analyze and visualize Hi-C data, Python. https://hicexplorer.readthedocs.io/en/latest/, https://github.com/deeptools/HiCExplorer/

- `hiclib` - Python tools to qc, map, normalize, filter and analyze Hi-C data, https://bitbucket.org/mirnylab/hiclib

- `HiC-Pro` - Python and command line-based optimized and flexible pipeline for Hi-C data processing, https://github.com/nservant/HiC-Pro
    - Servant, Nicolas, Nelle Varoquaux, Bryan R. Lajoie, Eric Viara, Chong-Jian Chen, Jean-Philippe Vert, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiC-Pro: An Optimized and Flexible Pipeline for Hi-C Data Processing.” Genome Biology 16 (December 1, 2015): 259. https://doi.org/10.1186/s13059-015-0831-x. - HiC pipeline, references to other pipelines, comparison. From raw reads to normalized matrices. Normalization methods, fast and memory-efficient implementation of iterative correction normalization (ICE). Data format. Using genotyping information to phase contact maps.

- `HiCUP` - Perl-based pipeline, alignment only, output - BAM files. http://www.bioinformatics.babraham.ac.uk/projects/hicup/
    - Wingett, Steven, Philip Ewels, Mayra Furlan-Magaril, Takashi Nagano, Stefan Schoenfelder, Peter Fraser, and Simon Andrews. “HiCUP: Pipeline for Mapping and Processing Hi-C Data.” F1000Research 4 (2015): 1310. https://doi.org/10.12688/f1000research.7334.1. - HiCUP pipeline, alignment only, removes artifacts (religations, duplicate reads) creating BAM files. Details about Hi-C sequencing artefacts. Used in conjunction with other pipelines.

- `HiTC` - R package for High Throughput Chromosome Conformation Capture analysis, https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. https://bioconductor.org/packages/release/bioc/html/HiTC.html and https://www.bioconductor.org/packages/devel/bioc/vignettes/HiTC/inst/doc/HiTC.pdf 

- `my5C`- web-based tools, well-documented analysis and visualization of 5S data, http://my5c.umassmed.edu/

