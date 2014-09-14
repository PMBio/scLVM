import h5py
import sys
#LIMIX
sys.path.append('../')
sys.path.append('./')
sys.path.append('../CFG')
sys.path.append('../include')
import scipy as SP
import pdb
import matplotlib as mpl
mpl.use('Agg')
import pylab as PL
import glob
import re
from utils import *
from barplot import *
from default import *

if __name__ == '__main__':
	from tcell import *
 
 #1. load data
	correlations=1#collect data from correlation analysis
	f = h5py.File(CFG['data_file'],'r')
	fpa = h5py.File(CFG['panama_file'],'r')
	cc_noise_filtered = fpa['cc_noise_filtered'][:]
	
	
	K = fpa['Kconf'][:]
	genes_het=SP.array(SP.where(f['genes_heterogen'][:].ravel()==1)).ravel()
	Nhet = len(genes_het)
	Ncells = f['LogNcountsMmus'][:].shape[0]
	Ngenes= f['LogNcountsMmus'][:].shape[1]
	
	cc_genes=SP.ones((Ngenes,1))
	cc_genes[cc_noise_filtered]=0
	cc_genes_het=cc_genes[genes_het]
	
	out_name = 'correlation_test'
	filename = 'Tcells'
	out_dir  = os.path.join(CFG['out_base'],out_name)
	run_dir  = os.path.join(out_dir,'runs'+filename)


	FL = glob.glob(os.path.join(run_dir,'job*.hdf5'))
	
	out_file_base = './vars'+filename
	out_file = out_file_base+'.hdf5' 
	fout = h5py.File(out_file,'w')

	RV = {}
	if correlations ==1:
		fields = ['pv','pv0','varsnorm','Ycorr','is_converged']
	else:
		fields = ['varsnorm','Ycorr','is_converged']

	for field in fields:
		RV[field] = SP.zeros([Nhet,Nhet])
	RV['varsnorm'] = SP.zeros([Nhet,4])
	RV['Ycorr'] = SP.zeros([Nhet,Ncells])
	RV['is_converged'] = SP.zeros([Nhet,1])

	#loop through files
	for fn in FL:
		#read file
		try:
			ff = h5py.File(fn,'r')
			for key in ff.keys():
				id0 = int(key.split('_')[1])
				for field in fields:
					if (field == 'Ycorr') or (field == 'varsnorm'):
						RV[field][id0,:] = ff[key][field][:]
					else:
						if (field == 'is_converged'):
							RV[field][id0,0] = ff[key][field][0,]*1.0
						else:
							RV[field][id0,:] = ff[key][field][:]
				pass
		except Exception:
			pdb.set_trace()
			print "bad file: %s" % (fn)
			continue
		
	#store
	#pdb.set_trace()
	RV['K'] = K
	RV['cc_genes_het'] = cc_genes_het
	dumpDictHdf5(RV,fout)
	fout.close()

	#plot variance decomposition			   
	from default import *
	indsconv = RV['is_converged'].ravel()==1 
	is_nocc_conv = SP.bitwise_and(indsconv,cc_genes_het.ravel()==1)
	vars_normConv=RV['varsnorm'][indsconv==1,:]
	vars_normConv = vars_normConv[:,[0,2,3,1]]
	H2=1-vars_normConv[:,2]
	out_file_pdf = out_file_base+'.pdf'
	var_plot(vars_normConv,H2,CFG['var_comp_fields'][SP.array([0,2,3,1])],filename=out_file_pdf,normalize=True, figsize=[5,4])	
