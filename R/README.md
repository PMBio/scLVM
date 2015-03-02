## scLVM and R

This folder contains 

* R scripts to process the filtered count table, both when spike-ins are present and for data-sets without spike-ins
  * A demo script which can be run for data where spike-ins are present: ``transform_counts_demo.Rmd``
  * A demo script showing how to estimate baseline variability, which is required by scLVM, for dat awithout spike-ins: ``transform_counts_demo_no_spikeins.Rmd``
*  Filtered count tables
  * ``data_Tcells.Rdata`` contains the T-cell data [1]
  * ``data_mESCquartz.Rdata`` contains the FPMK normalized mESC data (Quartz-Seq protocol) from  [2]
* An R wrapper wich is based on rPython.
  * A commented R markdown script illustrating how to use the R wrapper: ``scLVMr.Rmd``
  * The R wrapper ``R2py.R``
  * python helper functions ``init_data.py``

[1] Mahata, B., Zhang, X., Kolodziejczyk, A. A., Proserpio, V., Haim-Vilmovsky, L., Taylor, A. E., ... & Teichmann, S. A. (2014).
Single-cell RNA sequencing reveals T helper cells synthesizing steroids de novo to contribute to immune homeostasis. Cell reports, 7(4), 1130-1142.

[2] Sasagawa, Y., Nikaido, I., Hayashi, T., Danno, H., Uno, K. D., Imai, T., & Ueda, H. R. (2013). Quartz-Seq: a highly reproducible and sensitive single-cell RNA sequencing method, reveals non-genetic gene-expression heterogeneity. Genome Biol, 14(4), R31.

