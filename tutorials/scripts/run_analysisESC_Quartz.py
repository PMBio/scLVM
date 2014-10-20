import h5py
import sys
sys.path.append('./')
sys.path.append('../CFG')
sys.path.append('../include')
#limix_path = '/Users/florian/Code/python_code/limix-master/build/release.darwin/interfaces/python'
#sys.path.append(limix_path)
sys.path.append('./..')
#sys.path.append('../scLVM')
from scLVM import scLVM
import limix.modules.panama as PANAMA
import limix.modules.varianceDecomposition as VAR
from include.utils import dumpDictHdf5
import scipy as SP
import limix
import limix.modules.qtl as QTL
import limix.modules.varianceDecomposition as VAR
import pdb
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pylab as PL




if __name__ == '__main__':
	from quartz import *

	#Either distribute over several jobs 
	#(recommended when performing correlation analyis) or not (debug mode)
	if 'debug' in sys.argv:
		Njobs = 1
		ijob  = 0
	else:
		Njobs = int(sys.argv[1])
		ijob  = int(sys.argv[2])-1
	
	#Where to save results	
	out_name = 'correlation_test'
	out_dir  = os.path.join(CFG['out_base'],out_name)
	run_dir  = os.path.join(out_dir,'runsESCs')

	if not os.path.exists(run_dir):
		os.makedirs(run_dir)

	#load data
	f = h5py.File(CFG['data_file'],'r')
	Y = f['LogNcountsQuartz'][:]
	tech_noise = f['LogVar_techQuartz_logfit'][:]
	genes_het_bool=f['genes_heterogen'][:]	 # index of heterogeneous(??!??) genes
	geneID = f['gene_names_all'][:]			# gene names
	cellcyclegenes_filter = SP.unique(f['ccGO_gene_indices'][:].ravel() -1) # idx of cell cycle genes
	cellcyclegenes_filterCB600 = f['ccCBall_gene_indices'][:].ravel() -1		# idxof cell cycle genes ...
   

	# filter cell cycle genes
	idx_cell_cycle = SP.union1d(cellcyclegenes_filter,cellcyclegenes_filterCB600)
	Ymean2 = Y.mean(0)**2>0
	idx_cell_cycle_noise_filtered = SP.intersect1d(idx_cell_cycle,SP.array(SP.where(Ymean2.ravel()>0)))
	Ycc = Y[:,idx_cell_cycle_noise_filtered]
	
	#Fit GPLVM to data 
	k = 1					 # number of latent factors
	file_name = CFG['panama_file']# name of the cache file
	recalc = True # recalculate X and Kconf
	sclvm = scLVM(Y)
	X,Kcc,varGPLVM = sclvm.fitGPLVM(idx=idx_cell_cycle_noise_filtered,k=1,out_dir='./cache',file_name=file_name,recalc=recalc)

	#3. load relevant dataset for analysis
	genes_het=SP.array(SP.where(f['genes_heterogen'][:].ravel()==1))

   # considers only heterogeneous genes
	Ihet = genes_het_bool==1
	Y	= Y[:,Ihet]
	tech_noise = tech_noise[Ihet]
	geneID = geneID[Ihet] 
	
	
	#4. split across genes
	Iy	= SP.array(SP.linspace(0,Y.shape[1],Njobs+1),dtype='int')
	i0	= Iy[ijob]
	i1	= Iy[ijob+1]

	#create outfile
	out_file = os.path.join(run_dir,'job_%03d_%03d.hdf5' % (ijob,Njobs))
	fout	 = h5py.File(out_file,'w')
	
	#ground truth
	KS = f['KS'][:]
	KG1 = f['KG1'][:]
	KG2M = f['KG2M'][:]
	#KList = {}
	#KList[0] = KG1
	#KList[1] = KS
	#KList[2] = KG2M

	sclvm = scLVM(Y,geneID=geneID,tech_noise=tech_noise)

	# fit the model from i0 to i1
	sclvm.varianceDecomposition(K=Kcc,i0=i0,i1=i1) 
	normalize=True	# variance components are normalizaed to sum up to one

	# get variance components
	var, var_info = sclvm.getVarianceComponents(normalize=normalize)
#	var_filtered = var[var_info['conv']] # filter out genes for which vd has not converged


	# get corrected expression levels
	Ycorr = sclvm.getCorrectedExpression()
	
	# fit lmm without correction
	pv0,beta0,info0 = sclvm.fitLMM(K=None,i0=i0,i1=i1,verbose=True)
	# fit lmm with correction
	pv,beta,info = sclvm.fitLMM(K=Kcc,i0=i0,i1=i1,verbose=True)
	
	#write to file
	count = 0
	for i in xrange(i0,i1):
		gene_id = 'gene_%d' % (i)
		out_group = fout.create_group(gene_id)
		RV = {}
		RV['pv0'] = pv0[count,:]
		RV['pv'] = pv[count,:]
		RV['beta'] = beta[count,:]
		RV['beta0'] = beta0[count,:]
		RV['vars'] = var[count,:]
		RV['varsnorm'] = var[count,:]
		RV['is_converged']=SP.repeat(var_info['conv'][count]*1,5)
		RV['Ycorr'] = Ycorr[:,count]
		dumpDictHdf5(RV,out_group)
		count+=1
	fout.close()
		 
	
