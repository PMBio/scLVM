scLVM
What is scLVM?

scLVM is a modelling framework for single-cell RNA-seq data that can be used to dissect the observed heterogeneity into different sources, thereby allowing for the correction of confounding sources of variation. 

By Florian Buettner, Paolo Casale and Oliver Stegle (and ?)

Philosophy

Observed heterogeneity in single-cell profiling data is multi-factorial. scLVM provides an efficient framework for unravelling this heterogeneity, correcting for confounding facotrs and facilitating unbiased downstream analyses. scLVM builds on the Gaussian process latent variable mixed linear models. Our modelling approach is based on efficient inference algorithms implemented in [limix](https://github.com/PMBio/limix).

Installation:

* scLVM is particularly easy to install using the anaconda [anaconda](https://store.continuum.io/cshop/anaconda) python distribution.
 
* It requires Python 2.7 with
  - scipy, h5py, numpy, pylab

* scLVM relies heavily on [limix](https://github.com/PMBio/limix), which can be installed using ``pip install limix`` on most systems.

* If you would like to use the non-linear GPLVM for visualisation, you require the GPy package. This can be installed using `pip install GPy` 

* Preprocessing setps are executed in R and require R>3.0:

How to use scLVM?

A good starting point is the ipython notebook `tcells_demo.ipynb` in the ./demo folder. You can view a html version of the notebook or open it `using ipython notebook` and reproduce the results from Buettner et al 2014 [2]. The notebook requires an hdf5 data file containing the relevant data (normalised read counts etc.). We generate this data structure in R, and illustrate in `transform_counts_demo.Rmd` how this was done for the T-cell data. 


Problems ?

If you want to use LIMIX and encounter any issues, please contact us by email: fbuettner.phys@gmail.com

License

See 
