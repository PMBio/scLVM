"""configuration file"""
"""requires: base_dir, data_file and panama_file
base_dir   : where to find panama file and data file
data_file  : hdf5 file with the necessary data:

LogNcountsMmus       : log10(normalised reads)
cellcyclegenes_filter: indices of cell cycle genes, starting from 1 (eg from GO or cyclebase; filtered for expressed genes)
ccCBall_gene_indices : indices from separate source, eg cyclebase [optional]
LogVar_techMmus      : technical variance for each gene in log-space (eg estimated via spike-ins)
genes_heterogen      : boolean vector indicating whether a gene is variable above technical noise (e.g. as estimted following Brennecke et al)
panama_file: name of hdf5 file where cell-cycle induced convariance matrix is stored (generated with fit_panama.py)
"""
import os

CFG = {}
CFG['base_dir'] = '/Users/florian/Code/python_code/scLVM-1.0'
CFG['out_base'] = os.path.join(CFG['base_dir'],'out')
CFG['data_file'] = os.path.join(CFG['base_dir'],'data/mESC_C1/normCountsESC_Filter.h5f')
CFG['panama_file2d'] = os.path.join(CFG['out_base'],'Kkedarcc2d.h5f')
CFG['panama_file'] = os.path.join(CFG['out_base'],'panama_ESCcc.hdf5')
CFG['summary_file'] = os.path.join(CFG['out_base'],'varsESC.hdf5')
CFG['block_file'] = os.path.join(CFG['out_base'],'varsESCBlock.hdf5')
