#Tutorials


##scLVM
scLVM is implemented in python and can either be run it in python as demonstrated in the [demo ipython notebook](http://nbviewer.ipython.org/github/pmbio/scLVM/blob/master/tutorials/tcell_demo.ipynb) or as an R package as demonstrated [here](https://github.com/PMBio/scLVM/blob/master/R/tutorials/scLVM_vignette.Rmd).


In case you would like to distribute the computations over many cores, we provide scripts for easy parallelisation here: [scripts](https://github.com/PMBio/scLVM/tree/master/tutorials/scripts)


##Input
Our software requires as input a count table with read counts for all cells. This can be generated e.g. with [HTSeq](http://www-huber.embl.de/HTSeq/doc/overview.html) and appropriate filtering and low-level processing should be perfromed to filter out spurious/'bad' cells before the data are analyzed within scLVM.

Our [R package](https://github.com/PMBio/scLVM/tree/master/R) has convenience functions for processing this count table and the subsequent analysis with scLVM.

If you would like to use the python implementation, you can use our R scripts to process this filtered count table, both when spike-ins are present and for data-sets without spike-ins:

A demo script which can be run for data where spike-ins are present can be found here: [transform_counts_demo.Rmd](https://github.com/PMBio/scLVM/blob/master/R/scripts/transform_counts_demo.Rmd)

Without spike-ins, baseline variability, which is required by scLVM, can be estimated as shown here: [transform_counts_demo_no_spikeins.Rmd](https://github.com/PMBio/scLVM/blob/master/R/scripts/transform_counts_demo_no_spikeins.Rmd)

These R scripts generate an hdf5 file containng the variables required by the core scLVM algorithm:

* Normalised gene expression data -
* Technical noise (in log space)
* Gene symbols
* Heterogeneous genes (boolean vector)
* Cell cycle genes (vector of indices)
