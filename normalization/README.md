# Normalization of Hi-C data

- `HiCNorm` - removing biases in Hi-C data via Poisson regression, http://www.people.fas.harvard.edu/~junliu/HiCNorm/

- `HiTC` - The HiTC R package was developed to explore high-throughput 'C' data such as 5C or Hi-C. Dedicated R classes as well as standard methods for quality controls, normalization, visualization, and further analysis are also provided. https://bioconductor.org/packages/release/bioc/html/HiTC.html
    - Servant, Nicolas, Bryan R. Lajoie, Elphège P. Nora, Luca Giorgetti, Chong-Jian Chen, Edith Heard, Job Dekker, and Emmanuel Barillot. “HiTC: Exploration of High-Throughput ‘C’ Experiments.” Bioinformatics (Oxford, England) 28, no. 21 (November 1, 2012): 2843–44. https://doi.org/10.1093/bioinformatics/bts521. - HiTC paper. Processed data import from TXT/BED into GRanges. Quality control, visualization. Normalization using loess regression on genomic distance, 45-degree rotation and visualization of triangle TADs. Adding annotation at the bottom. PCA to detect A/B compartments. 

- `caICB` - chromosome-adjusted iterative correction of biases method. Review of Hi-C data biases and methods to remove them. CNV bias cannot be removed after within-chromosome iterative correction bias methods. The distance-dependent decay of interaction frequencies is modeled by splines. The caICB method aims to minimize the differences across count-distance curves of different chromosomes, chr1 as the reference. A priori knowledge of CNVs is not required. The between- and within-chromosome bias is removed, minimizing the number of significant contacts due to CNVs. https://bitbucket.org/mthjwu/hicapp
    - Wu, Hua-Jun, and Franziska Michor. “A Computational Strategy to Adjust for Copy Number in Tumor Hi-C Data.” Bioinformatics (Oxford, England) 32, no. 24 (15 2016): 3695–3701. https://doi.org/10.1093/bioinformatics/btw540.

- `ICE` - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. 
    - Imakaev, Maxim, Geoffrey Fudenberg, Rachel Patton McCord, Natalia Naumova, Anton Goloborodko, Bryan R. Lajoie, Job Dekker, and Leonid A. Mirny. “Iterative Correction of Hi-C Data Reveals Hallmarks of Chromosome Organization.” Nature Methods 9, no. 10 (October 2012): 999–1003. https://doi.org/10.1038/nmeth.2148. - ICE - Iterative Correction and Eigenvalue decomposition, normalization of HiC data. Assumption - all loci should have equal visibility. Deconvolution into eigenvectors/values. hiclib https://bitbucket.org/mirnylab/hiclib. Good description of the algorithm by Lior Pachter https://liorpachter.wordpress.com/2013/11/17/imakaev_explained/

- `OneD` - CNV bias-correction method, addresses the problem of partial aneuploidy. Bin-centric counts are modeled using negative binomial distribution, and its parameters are estimated using splines. A hidden Markov model is fit to infer copy number for each bin. Each Hi-C matrix entry is corrected by dividing its value by square root of the product of CNVs for the corresponding bins. Reproducibility score (eigenvector decomposition and comparison) to measure improvement in the similarity between replicated Hi-C data. https://github.com/qenvio/dryhic
    - Vidal, Enrique, François le Dily, Javier Quilez, Ralph Stadhouders, Yasmina Cuartero, Thomas Graf, Marc A Marti-Renom, Miguel Beato, and Guillaume J Filion. “OneD: Increasing Reproducibility of Hi-C Samples with Abnormal Karyotypes.” Nucleic Acids Research, January 31, 2018. https://doi.org/10.1093/nar/gky064.



## Papers

- Yaffe, Eitan, and Amos Tanay. “Probabilistic Modeling of Hi-C Contact Maps Eliminates Systematic Biases to Characterize Global Chromosomal Architecture.” Nature Genetics 43, no. 11 (November 2011): 1059–65. https://doi.org/10.1038/ng.947. - Sources of biases: 1) non-specific ligation (large distance between pairs); 2) length of each ligated fragments; 3) CG content and nucleotide composition; 4) Mappability. Normalization. Enrichment of long-range interactions in active promoters. General aggregation of active chromosomal domains. Chromosomal territories, high-activity and two low-activity genomic clusters
