# A collection of tools using Hi-C data 

## Reviews

- A list of tools available for data analysis and/or visualization of 4DN-related datasets. https://www.4dnucleome.org/software.html

- Ay, Ferhat, and William S. Noble. “Analysis Methods for Studying the 3D Architecture of the Genome.” Genome Biology 16 (September 2, 2015): 183. https://doi.org/10.1186/s13059-015-0745-7. - Hi-C technology and methods review. [Table 1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0745-7#Tab1) - list of tools. Biases, normalization, matrix balancing. Extracting significant contacts, obs/exp ratio, parametric (powerlaw, neg binomial, double exponential), non-parametric (splines). 3D enrichment. References. TAD identification, directionality index. Outlook, importance of comparative analysis

- Forcato, Mattia, Chiara Nicoletti, Koustav Pal, Carmen Maria Livi, Francesco Ferrari, and Silvio Bicciato. “Comparison of Computational Methods for Hi-C Data Analysis.” Nature Methods, June 12, 2017. https://doi.org/10.1038/nmeth.4325. - Hi-C processing and TAD calling tools benchmarking, [Table 1](https://www.nature.com/articles/nmeth.4325/tables/1), simulated (Lun and Smyth method) and real data. Notes about pluses and minuses of each tool. TAD reproducibility is higher than chromatin interactions, increases with larger number of reads. Consistent enrichment of TAD boundaries in CTCF, irrespectively of TAD caller. Hi-C replication is poor, just a bit more than random. Supplementary table 2 - technical details about each program, Supplementary Note 1 - Hi-C preprocessing tools, Supplementary Note 2 - TAD callers. Supplementary note 3 - how to simulate Hi-C data. Supplementary note 6 - how to install tools. https://images.nature.com/full/nature-assets/nmeth/journal/v14/n7/extref/nmeth.4325-S1.pdf

- Nicoletti, Chiara, Mattia Forcato, and Silvio Bicciato. “Computational Methods for Analyzing Genome-Wide Chromosome Conformation Capture Data.” Current Opinion in Biotechnology 54 (December 2018): 98–105. https://doi.org/10.1016/j.copbio.2018.01.023. - 3C-Hi-C tools review, Table 1 lists categorizes main tools, Figure 1 displays all steps in technology and analysis (alignment, resolution, normalization, including accounting for CNVs, A/B compartments, TAD detection, visualization). Concise description of all tools.



## Multi-purpose

- `HiCExplorer` - set of programs to process, normalize, analyze and visualize Hi-C data, Python. https://hicexplorer.readthedocs.io/en/latest/, https://github.com/deeptools/HiCExplorer/

- `hiclib` - tools to qc, map, normalize, filter and analyze Hi-C data, https://bitbucket.org/mirnylab/hiclib

- `my5C` tools - well-documented analysis and visualization of 5S data, http://my5c.umassmed.edu/


## Visualization

- `HiTC` - import, normalize, annotate and visualize Hi-C data, https://bioconductor.org/packages/release/bioc/html/HiTC.html and https://www.bioconductor.org/packages/devel/bioc/vignettes/HiTC/inst/doc/HiTC.pdf

- `pyGenomeTracks` - Standalone program and library to plot beautiful genome browser tracks. https://github.com/deeptools/pyGenomeTracks


## Normalization

- `HiCNorm` - removing biases in Hi-C data via Poisson regression, http://www.people.fas.harvard.edu/~junliu/HiCNorm/


## Interaction callers

- `Fit-Hi-C` - detection of significant chromatin interactions, https://noble.gs.washington.edu/proj/fit-hi-c/


## TAD detection

- `Armatus` - TAD detection at different resolutions, https://www.cs.cmu.edu/~ckingsf/software/armatus/, https://github.com/kingsfordgroup/armatus

- `HiCseg` - TAD detection by maximization of likelihood based block-wise segmentation model, https://cran.r-project.org/web/packages/HiCseg/index.html

## QC

- `3DChromatin_ReplicateQC` - Yardimci, Galip, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data,” September 14, 2017. doi:10.1101/188755. - Comparison of four Hi-C reproducibility assessment tools, `HiCRep`, `GenomeDISCO`, `HiC-Spector`, `QuASAR-Rep`. Tested the effects of noise, sparsity, resolution. Spearman doesn't work well. All tools performed similarly, worsening expectedly. QuASAR has QC tool measuring the level of noise. https://github.com/kundajelab/3DChromatin_ReplicateQC


## SNP-oriented

- `iRegNet3D` - Integrated Regulatory Network 3D (iRegNet3D) is a high-resolution regulatory network comprised of interfaces of all known transcription factor (TF)-TF, TF-DNA interaction interfaces, as well as chromatin-chromatin interactions and topologically associating domain (TAD) information from different cell lines.  
    - Goal: SNP interpretation
    - Input: One or several SNPs, rsIDs or genomic coordinates.
    - Output: For one or two SNPs, on-screen information of their disease-related info, connection over TF-TF and chromatin interaction networks, and whether they interact in 3D and located within TADs. For multiple SNPs, same info downloadable as text files.

