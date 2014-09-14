"""configuration file"""
"""requires: base_dir, data_file and panama_file
    base_dir   : where to find panama file and data file
    data_file  : hdf5 file with the necessary data:
    LogNcountsMmus       : log10(normalised reads)
    cellcyclegenes_filter: indices of cell cycle genes, starting from 1 (eg from GO or cyclebase; filtered for expressed genes)
    LogVar_techMmus      : technical variance for each gene in log-space (eg estimated via spike-ins)
    genes_heterogen      : boolean vector indicating whether a gene is variable above technical noise (e.g. as estimted following Brennecke et al)
    panama_file: name of hdf5 file where cell-cycle induced convariance matrix is stored (generated with fit_panama.py)
"""

import os


CFG = {}
CFG['base_dir'] = '/Users/florian/Code/limixFlorian/python/scLVM'
CFG['out_base'] = os.path.join(CFG['base_dir'],'out')
CFG['summary_file'] = os.path.join(CFG['out_base'],'correlation_test/runsQuartz3Klmm.hdf5')
CFG['data_file'] = os.path.join(CFG['base_dir'],'data/quartz/normCounts_mESCquartz.h5f')
CFG['panama_file'] = os.path.join(CFG['out_base'],'panama_quartz_GOCB.hdf5')
