# The Smooth Euler Characteristic Transform

Radiomics is focused on the extraction of quantitative features from medical images, typically constructed by tomography and digitally stored as shapes or surfaces. In [Crawford et al. (2017)](https://arxiv.org/abs/1611.06818), we explore the use of a novel statistic, the smooth Euler characteristic transform (SECT), as an automated procedure to extract geometric or topological statistics from tumor images. More generally, the SECT is designed to integrate shape information into regression models by representing shapes and surfaces as a collection of curves. Due to its well-defined inner product structure, the SECT can be used in a wider range of functional and nonparametric modeling approaches than other previously proposed topological summary statistics. We illustrate the utility of the SECT in a radiomics context by showing that the topological quantification of tumors, assayed by magnetic resonance imaging (MRI), are better predictors of clinical outcomes in patients with glioblastoma multiforme (GBM). Using publicly available data from [The Cancer Genome Atlas (TCGA)](https://portal.gdc.cancer.gov) and [The Cancer Imaging Archive (TCIA)](https://wiki.cancerimagingarchive.net/display/Public/TCGA-GBM), we show that SECT features alone explain more of the variance in patient survival than gene expression, volumetric features, and morphometric features.

The smooth Euler characteristic transform for images and three-dimensional shapes is implemented as a set of R and MATLAB routines. Results in the manuscript were derived by using the R version of the code. The Bayesian Gaussian process (GP) regression model that we used to incorporate shape statistics in predictive analyses is also carried out within the R environment. 

### The MATLAB Environment
MATLAB is a multi-paradigm numerical computing environment and a proprietary programming language developed by [MathWorks](https://www.mathworks.com/index-c.html). MATLAB allows matrix manipulations, plotting of functions and data, implementation of algorithms, creation of user interfaces, and interfacing with programs written in other languages. For more on licensing options, please visit [here](https://www.mathworks.com/campaigns/products/ppc/google/matlab-toolbox-price-request.html?form_seq=reg).

### The R Environment
R is a widely used, free, and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the [Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/) CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions that are likely sufficient for many users' needs. Users can also install R from source code;  however, this may require a significant amount of effort. For specific details on how to compile, install, and manage R and R-packages, refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).

### R Packages Required for SECT and the GP Regression
The statistical implementation of the SECT topological summaries using functional RKHS regression models requires the installation of the following R libraries:

[BGLR](https://cran.r-project.org/web/packages/BGLR/index.html)

[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)

[gdata](https://cran.r-project.org/web/packages/gdata/index.html)

[Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)

[RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/index.html)

[RcppParallel](https://cran.r-project.org/web/packages/RcppParallel/index.html)

[BAKR](https://github.com/lorinanthony/BAKR)

[R.matlab](https://cran.r-project.org/web/packages/R.matlab/index.html) (If using the MATLAB version of SECT)

The easiest method to install these packages is with the following example command entered in an R shell:

    install.packages("BGLR", dependecies = TRUE)

Alternatively, one can also [install R packages from the command line](http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages).

### C++ Functions Required for RKHS Software Packages
The code in this repository assumes that basic C++ functions and applications are already set up on the running personal computer or cluster. If not, some of the BAKR/BGLR functions and necessary Rcpp packages will not work properly. A simple option is to use [gcc](https://gcc.gnu.org/). macOS users may use this collection by installing the [Homebrew package manager](http://brew.sh/index.html) and then typing the following into the terminal:

    brew install gcc

For macOS users, the Xcode Command Line Tools include a GCC compiler. Instructions on how to install Xcode may be found [here](http://railsapps.github.io/xcode-command-line-tools.html). For extra tips on how to run C++ on macOS, please visit [here](http://seananderson.ca/2013/11/18/rcpp-mavericks.html). For tips on how to avoid errors dealing with "-lgfortran" or "-lquadmath", please visit [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

### Data Availability
The results shown here are in whole or part based upon data generated by the TCGA Research Network. DICOM formatted MRI scans and patient clinical information were taken directly from the TCIA web portal. Matched molecular data were downloaded directly from the Genomic Data Commons (GDC) by selecting the RNA-Seq tab option under [GBM](https://portal.gdc.cancer.gov/projects/TCGA-GBM). Shape-based summary statistics necessary for replicating this study (i.e. the segmented tumor images, the volumetric measurements, morphometric data, and topological summary statistics) are publicly available on this repository under Data.

### Segmented TCIA Magnetic Resonance Images (MRIs)
MRIs of primary GBM tumors were collected from 90+ patients archived by the TCIA, which is a publicly accessible data repository of medical images of cancer patients with matched data in TCGA. Briefly, TCGA provides a collection of a variety of genomic and clinical data for 33 types of cancer. The patients analyzed in Crawford et al. (2017) were selected based on two sets of criteria, namely: (1) they had post-contrast T1 axial MRIs taken at the time of their diagnosis, and (2) they had available matching (mRNA) gene expression data and clinical correlates on [cBioPortal](http://www.cbioportal.org).

We segmented the TCIA MRIs using the [Medical Imaging Interaction Toolkit with augmented tools for segmentation (MITKats)](https://github.com/RabadanLab/MITKats), from [Chen and Rabadan (2017)](https://link.springer.com/chapter/10.1007/978-3-319-69775-8_10), to extract tumor lesions from the surrounding brain tissue. Briefly, this algorithm first converts MRI images to a grayscale, and then thresholds to generate binary images. Morphological segmentation is then applied to delineate connected components. More specifically, the program selects contours corresponding to enhanced tumor lesions, which are lighter than healthy brain tissue. For instance, necrosis presents as dark regions nested within the indicated lesion. Examples of the raw image obtained from TCIA, and the final segmented result, are given in Crawford et al. (2017) under Figure 3(a) and Figure 3(b), respectively. All segmented TCIA images used in our study can be found in a zipped file in the Data directory.

### Running SECT and GP Regression
The tutorial for running a predictive analysis, similar to the one presented in Crawford et al. (2017), can be found in the Analysis directory. Note that the current version of the SECT code takes images/shapes that are formatted as png files. This code looks at a subset of the TCGA patients. This script serves as a means of reproducibility of the results presented in the manuscript.

### Relevant Citations
L. Crawford, A. Monod, A.X. Chen, S. Mukherjee, and R. Rabadán (2017). Functional Data Analysis using a Topological Summary Statistic: the Smooth Euler Characteristic Transform. arXiv. 1611.06818.

### Questions and Feedback
For questions or concerns with the SECT functions, please contact [Lorin Crawford](mailto:lorin_crawford@brown.edu) or [Anthea Monod](mailto:am4691@cumc.columbia.edu).

We appreciate any feedback you may have with our repository and instructions.
