## Reproducibility and QC of Hi-C data

- `3DChromatin_ReplicateQC` - Comparison of four Hi-C reproducibility assessment tools, `HiCRep`, `GenomeDISCO`, `HiC-Spector`, `QuASAR-Rep`. Tested the effects of noise, sparsity, resolution. Spearman doesn't work well. All tools performed similarly, worsening expectedly. QuASAR has QC tool measuring the level of noise. https://github.com/kundajelab/3DChromatin_ReplicateQC
    - Yardimci, Galip, Hakan Ozadam, Michael E.G. Sauria, Oana Ursu, Koon-Kiu Yan, Tao Yang, Abhijit Chakraborty, et al. “Measuring the Reproducibility and Quality of Hi-C Data,” September 14, 2017. doi:10.1101/188755. 

- `HiC-Spector` - reproducibility metric to quantify the similarity between contact maps using spectral decomposition. Decomposing Laplacian matrices and sum the Euclidean distance between eigenvectors. https://github.com/gersteinlab/HiC-spector
    - Yan, Koon-Kiu, Galip Gürkan Yardimci, Chengfei Yan, William S. Noble, and Mark Gerstein. “HiC-Spector: A Matrix Library for Spectral and Reproducibility Analysis of Hi-C Contact Maps.” Bioinformatics (Oxford, England) 33, no. 14 (July 15, 2017): 2199–2201. https://doi.org/10.1093/bioinformatics/btx152.

