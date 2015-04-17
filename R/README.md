## scLVM and R

This folder contains 

* An R package wich is based on rPython.
  * To install the package download the .tar.gz file and run R CMD INSTALL sclvm_0.99.1.tar.gz
  * Source code can be found in the scLVM folder
* R markdown scripts illustrating how to use the R interface for data with and without spike-ins. We also illustrate how to use scLVM in the presence of multiple facotrs. All scripts are located in the tutorial folder.
  * A demo script which can be run for data where spike-ins are present: ``scLVMr_demo.Rmd``
  *  Filtered count tables needed to run the exanmples are part of the package and located in the data folder
    * ``data_Tcells.rda`` contains the T-cell data [1]
    * ``data_mESCquartz.rda`` contains the FPMK normalized mESC data (Quartz-Seq protocol) from  [2]




[1] Mahata, B., Zhang, X., Kolodziejczyk, A. A., Proserpio, V., Haim-Vilmovsky, L., Taylor, A. E., ... & Teichmann, S. A. (2014).
Single-cell RNA sequencing reveals T helper cells synthesizing steroids de novo to contribute to immune homeostasis. Cell reports, 7(4), 1130-1142.

[2] Sasagawa, Y., Nikaido, I., Hayashi, T., Danno, H., Uno, K. D., Imai, T., & Ueda, H. R. (2013). Quartz-Seq: a highly reproducible and sensitive single-cell RNA sequencing method, reveals non-genetic gene-expression heterogeneity. Genome Biol, 14(4), R31.

