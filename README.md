# A collection of tools using Hi-C data 

## General

- `hiclib` - tools to qc, map, normalize, filter and analyze Hi-C data, https://bitbucket.org/mirnylab/hiclib

- `my5C` tools - well-documented analysis and visualization of 5S data, http://my5c.umassmed.edu/


## Visualization

- `HiTC` - import, normalize, annotate and visualize Hi-C data, https://bioconductor.org/packages/release/bioc/html/HiTC.html and https://www.bioconductor.org/packages/devel/bioc/vignettes/HiTC/inst/doc/HiTC.pdf

- `pyGenomeTracks` - Standalone program and library to plot beautiful genome browser tracks. https://github.com/deeptools/pyGenomeTracks


## Normalization

- `HiCNorm` - removing biases in Hi-C data via Poisson regression, http://www.people.fas.harvard.edu/~junliu/HiCNorm/


## QC

- `3DChromatin_ReplicateQC` - Yardimci, Galip, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data,” September 14, 2017. doi:10.1101/188755. - Comparison of four Hi-C reproducibility assessment tools, `HiCRep`, `GenomeDISCO`, `HiC-Spector`, `QuASAR-Rep`. Tested the effects of noise, sparsity, resolution. Spearman doesn't work well. All tools performed similarly, worsening expectedly. QuASAR has QC tool measuring the level of noise. https://github.com/kundajelab/3DChromatin_ReplicateQC


## SNP-oriented

- `iRegNet3D` - Integrated Regulatory Network 3D (iRegNet3D) is a high-resolution regulatory network comprised of interfaces of all known transcription factor (TF)-TF, TF-DNA interaction interfaces, as well as chromatin-chromatin interactions and topologically associating domain (TAD) information from different cell lines.  
    - Goal: SNP interpretation
    - Input: One or several SNPs, rsIDs or genomic coordinates.
    - Output: For one or two SNPs, on-screen information of their disease-related info, connection over TF-TF and chromatin interaction networks, and whether they interact in 3D and located within TADs. For multiple SNPs, same info downloadable as text files.

