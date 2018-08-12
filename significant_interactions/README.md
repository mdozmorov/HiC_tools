## Interaction callers

- `Fit-Hi-C` - Python tool for detection of significant chromatin interactions, https://noble.gs.washington.edu/proj/fit-hi-c/
    - Ay, Ferhat, Timothy L. Bailey, and William Stafford Noble. “Statistical Confidence Estimation for Hi-C Data Reveals Regulatory Chromatin Contacts.” Genome Research 24, no. 6 (June 2014): 999–1011. https://doi.org/10.1101/gr.160374.113. - Fit-Hi-C method, Splines to model distance dependence. Model mid-range interaction frequencies, decay with distance. Biases, methods for normalization. Two-step splines - use all dots for first fit, identify and remove outliers, second fit without outliers. Markers of boundaries - insulators, heterochromatin, pluripotent factors. CNVs are enriched in chromatin boundaries. Replication timing data how-to http://www.replicationdomain.com/. Validation Hi-C data. http://chromosome.sdsc.edu/mouse/hi-c/download.html

- `GoTHIC` - R package for peak calling in individual HiC datasets, while accounting for noise. https://www.bioconductor.org/packages/release/bioc/html/GOTHiC.html
    - Mifsud, Borbala, Inigo Martincorena, Elodie Darbo, Robert Sugar, Stefan Schoenfelder, Peter Fraser, and Nicholas M. Luscombe. “GOTHiC, a Probabilistic Model to Resolve Complex Biases and to Identify Real Interactions in Hi-C Data.” Edited by Mark Isalan. PLOS ONE 12, no. 4 (April 5, 2017): e0174744. https://doi.org/10.1371/journal.pone.0174744. - The GOTHiC (genome organisation through HiC) algorithm uses a simple binomial distribution model to simultaneously remove coveralge-associated biases in Hi-C data and detect significant interactions by assuming that the global background interaction frequency of two loci. Use of the Benjamini–Hochberg multiple-testing correction to control for the false discovery rate. 

- `HOMER` - Perl scripts for normalization, visualization, significant interaction detection, motif discovery. Does not correct for bias. http://homer.ucsd.edu/homer/interactions/

