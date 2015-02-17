#scLVM


##What is scLVM?

scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation. 

scLVM was primarily designed to account for cell-cycle induced variations in single-cell RNA-seq data where cell cycle is the primary soure of variability. For other use cases tutorials will follow shortly.

Software by Florian Buettner, Paolo Casale and Oliver Stegle. scLVM is explained in more detail in the accompanying publication [1].

##Philosophy

Observed heterogeneity in single-cell profiling data is multi-factorial. scLVM provides an efficient framework for unravelling this heterogeneity, correcting for confounding factors and facilitating unbiased downstream analyses. scLVM builds on Gaussian process latent variable models and linear mixed models. The underlying models are based on inference schemes implemented in [LIMIX](https://github.com/PMBio/limix).

##Installation:

* scLVM is particularly easy to install using the [anaconda](https://store.continuum.io/cshop/anaconda) python distribution. 
 
* It requires Python 2.7 with
  - scipy, h5py, numpy, pylab

* scLVM relies heavily on [limix](https://github.com/PMBio/limix) (version 0.6.4 or higher), which can be installed using ``pip install limix`` on most systems.

* If you would like to use the non-linear GPLVM for visualisation, we suggest installing the [GPy](https://github.com/SheffieldML/GPy) package. This can be installed using `pip install GPy`.

* Preprocessing steps are executed in R and require R>3.0:
For an example of how raw counts can be processed appropriately, see [R/transform_counts_Tcells.R](https://github.com/PMBio/scLVM/blob/master/R/transform_counts_Tcells.R). A more extensively commented markdown version of the script can be found here: [transform_counts_demo.Rmd](https://github.com/PMBio/scLVM/blob/master/R/transform_counts_demo.Rmd).

* For users who prefer to run the entire scLVM pipeline in R, we also provide an R wrapper wich is based on [rPython](http://cran.r-project.org/web/packages/rPython/index.html). A commented R markdown script can be found here: [scLVMr.Rmd](https://github.com/PMBio/scLVM/blob/master/R/scLVMr.Rmd).

##How to use scLVM?
The current software version should be considered as beta. Still, the method is working and can be used to reproduce the result of the accompanying publication [1]. More extensive documentation, tutorials and examples will be available soon.

A good starting point are the [tutorials](https://github.com/PMBio/scLVM/blob/master/tutorials).

For an illustration of how scLVM can be applied to the T-cell data considered in Buettner et al. [1], we have prepared a notebook that can be viewed [interactively](http://nbviewer.ipython.org/github/pmbio/scLVM/blob/master/tutorials/tcell_demo.ipynb) or alternatively as [PDF](https://github.com/PMBio/scLVM/blob/master/tutorials/tcell_demo.pdf) export.


##Problems ?

If you want to use scLVM and encounter any issues, please contact us by email: scLVM-dev@ebi.ac.uk

##License
See [LICENSE](https://github.com/PMBio/scLVM/blob/master/license.txt)

###References
[1] Buettner F, Natarajan KN, Casale FP, Proserpio V, Scialdone A, Theis FJ, Teichmann SA, Marioni JC & Stegle O, 2014. Computational analysis of cell-to-cell heterogeneity in single-cell RNA-Sequencing data reveals hidden subpopulation of cells, Nature Biotechnology, doi: 10.1038/nbt.3102.
