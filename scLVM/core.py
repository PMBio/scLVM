# Copyright(c) 2014, The scLVM developers (Forian Buettner, Paolo Francesco Casale, Oliver Stegle)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

import sys
sys.path.append('./..')
#import limix #make sure we use the right limix version

from utils.misc import dumpDictHdf5
from utils.misc import PCA 
from utils.misc import warning_on_one_line 
import limix
try:
    limix.__version__
    if versiontuple(limix.__version__)>versiontuple('0.7.3'):
        import limix.deprecated as limix

except:
    import limix.deprecated as limix    

import limix.modules.panama as PANAMA
import limix.modules.varianceDecomposition as VAR
import limix.modules.qtl as QTL

import scipy as SP
import scipy.linalg
import scipy.stats
import pdb
import h5py
import time
import copy
import warnings
import os

class scLVM:
	"""
	Single Cell Latent Varible Model module (scLVM)
	This class takes care of fitting and interpreting latent variable models to account for confounders in  single-cell RNA-Seq data 
	This module requires LIMIX
	"""

	def __init__(self,Y,geneID=None,tech_noise=None):
     
		"""
		Args:
			Y:			  	gene expression matrix [N, G]
			geneID:		  	G vector of geneIDs
			tech_noise:		G vector of tech_noise
		"""
		
		#store dimensions
		self.N = Y.shape[0]
		self.G = Y.shape[1]

		#set data
		self.Y	  = Y
		self.geneID = geneID
		self.tech_noise = None
		self.var=None
		if tech_noise is not None:
			self.set_tech_noise(tech_noise)

	
	def fitGPLVM(self,idx=None,k=1,standardize=False,out_dir='./cache',file_name=None,recalc=False, use_ard=False):
		"""
		Args:
			idx:			index of the genes involved
							(e.g., for capturing cell cycle, index of cell cycle genes)
			k:				number of latent factors
			standardize:	if True, rescale gene expression by std prior to fitting GPLVM
							(data are always mean-centered)
			out_dir:		dir used to cache the results
			file_name:		if not None, caches the results in the out_dir if the file does not exist
							if the file exists loads the results if recalc is True
			recalc:			if True and cache file exists, rewrite cacheFile
			use_ard:		use automatic relevance detection (switch off unimportant factors)
		Returns:
			X:				hidden variable
			Kconf:			similarity matrix based on the confounding effect (XX.T)
			varGPLVM:		variance contributions of latent factors and residual biological noise
		"""
		assert idx is not None, 'scLVM:: specify idx'
		if use_ard==True and  k<2:
			warnings.formatwarning = warning_on_one_line
			warnings.warn('when using ARD consider choosing k>1')

		file_out = os.path.join(out_dir,file_name)
		if not os.path.exists(file_out) or recalc==True:
			# prepare data
			Yconf = self.Y[:,idx]
			Yconf-= Yconf.mean(0)
			# fit gplvm
			panama = PANAMA.PANAMA(Y=Yconf,use_Kpop=False,standardize=standardize)
			panama.train(rank=k,LinearARD=use_ard)			
			X	 = panama.get_Xpanama()
			Kconf = panama.get_Kpanama()
			var = panama.get_varianceComps()
			if use_ard==False:
				varGPLVM = {'K':var['Kpanama'],'noise':var['noise']}
			else:
				varGPLVM = {'X_ARD':var['LinearARD'],'noise':var['noise']}
			# export results
			if not os.path.exists(out_dir):
				os.makedirs(out_dir)
			fout = h5py.File(file_out,'w')
			RV = {'X':X,'Kconf':Kconf}
			RV['cc_noise_filtered'] = idx
			dumpDictHdf5(RV,fout)
			dumpDictHdf5(varGPLVM,fout)
			fout.close()
		else:
			# load results from the file
			f = h5py.File(file_out,'r')
			X = f['X'][:]; Kconf = f['Kconf'][:]
			if use_ard==False:
				varGPLVM = {'K':f['K'][:],'noise':f['noise'][:]}
			else:	
				varGPLVM = {'X_ARD':f['X_ARD'][:],'noise':f['noise'][:]}
			f.close()

		return X,Kconf,varGPLVM

	def set_tech_noise(self,tech_noise):
		"""
		Args:
			tech_noise:		G vector of technical noise
		"""
		assert tech_noise.shape[0]==self.G, 'scLVM:: tech_noise dimension dismatch'
		
		self.tech_noise = tech_noise
		
	def varianceDecomposition(self,K=None,tech_noise=None,idx=None,i0=None,i1=None,max_iter=10,verbose=False):
		"""
		Args:
			K:				list of random effects to be considered in the analysis
			idx:			indices of the genes to be considered in the analysis
			i0:				gene index from which the anlysis starts
			i1:				gene index to which the analysis stops
			max_iter:		maximum number of random restarts
			verbose:		if True, print progresses
		"""

		if tech_noise is not None:		self.set_tech_noise(tech_noise)
		assert self.tech_noise is not None, 'scLVM:: specify technical noise'
		assert K is not None, 'scLVM:: specify K'

		if type(K)!=list:	K = [K]
		for k in K:
			assert k.shape[0]==self.N, 'scLVM:: K dimension dismatch'
			assert k.shape[1]==self.N, 'scLVM:: K dimension dismatch'

		if idx is None:
			if i0==None or i1==None:
				i0 = 0; i1 = self.G
			idx = SP.arange(i0,i1)
		elif type(idx)!=SP.ndarray:
			idx = SP.array([idx])

		_G	 = len(idx)
		var	= SP.zeros((_G,len(K)+2))
		_idx   = SP.zeros(_G)
		geneID = SP.zeros(_G,dtype=str)
		conv   = SP.zeros(_G)==1
		Ystar  = [SP.zeros((self.N,_G)) for i in range(len(K))]
		count  = 0
		Ystd = self.Y-self.Y.mean(0) #delta optimization might be more efficient
		Ystd/= self.Y.std(0)
		tech_noise = self.tech_noise/SP.array(self.Y.std(0))**2
		for ids in idx:
			if verbose:
				print '.. fitting gene %d'%ids
			# extract a single gene
			y = Ystd[:,ids:ids+1]
			# build and fit variance decomposition model
			vc= VAR.VarianceDecomposition(y)
			vc.addFixedEffect()
			for k in K:
				vc.addRandomEffect(k)
			vc.addRandomEffect(SP.eye(self.N))
			vc.addRandomEffect(SP.eye(self.N))
			vc.vd.getTerm(len(K)+1).getKcf().setParamMask(SP.zeros(1))
			for iter_i in range(max_iter):
				scales0 = y.std()*SP.randn(len(K)+2)
				scales0[len(K)+1]=SP.sqrt(tech_noise[ids]);
				_conv = vc.optimize(scales0=scales0)
				if _conv: break
			conv[count] = _conv
			if not _conv:
				var[count,-2] = SP.maximum(0,y.var()-tech_noise[ids])
				var[count,-1] = tech_noise[ids]
				count+=1;
				if self.geneID is not None:	geneID[count] = self.geneID[ids]
				continue
			_var = vc.getVarianceComps()[0,:]
			KiY = vc.gp.agetKEffInvYCache().ravel()
			for ki in range(len(K)):
				Ystar[ki][:,count]=_var[ki]*SP.dot(K[ki],KiY)
			var[count,:] = _var
			count+=1;
	
		# col header
		col_header = ['hidden_%d'%i for i in range(len(K))]
		col_header.append('biol_noise')
		col_header.append('tech_noise')
		col_header = SP.array(col_header)

		# annotate column and rows of var and Ystar
		var_info = {'gene_idx':idx,'col_header':col_header,'conv':conv}
		if geneID is not None:	var_info['geneID'] = SP.array(geneID)
		Ystar_info = {'gene_idx':idx,'conv':conv}
		if geneID is not None:	Ystar_info['geneID'] = SP.array(geneID)

		# cache stuff
		self.var   = var
		self.Ystar = Ystar
		self.var_info   = var_info
		self.Ystar_info = Ystar_info

	def getVarianceComponents(self,normalize=False):
		"""
		Returns:
			var:			variance component matrix [G_0,m]
							(m = number of variance components=len(K)+2)
							(G0 = genes which were considered)
			normalize:		if True, variance components are normalized to sum up to 1 (default False)
			var_info:		dictionary annotating rows and columns of var
							gene_idx:	index of the gene in the full matrix
							col_header: labels of variance components in the matrix
							conv:		boolean vetor marking genes for which variance decomposition has converged
							geneID:	 annotate rows of the variance component matrix
		"""
		assert self.var is not None, 'scLVM:: use varianceDecomposition method before'
		if normalize:	var = self.var/self.var.sum(1)[:,SP.newaxis]
		else:			var = self.var
		return var, self.var_info
 
	def getPredictions(self):
		"""
		Returns:
			Ystar:			predictions [N,G_0]
							(G0 = genes which were considered in the analysis)
							Remark: if more than 1 random effect is considered Ystar is a list of predicion matrices, ones per random effec
			Ystar_info:		annotate rows and columns of Ystar
							gene_idx:	index of the gene in the full matrix
							conv:		boolean vetor marking genes for which variance decomposition has converged
							geneID:	 annotate rows of the variance component matrix
		"""
		assert self.var is not None, 'scLVM:: use varianceDecomposition method before'
		if len(self.Ystar)==1:	Ystar = self.Ystar[0]
		else:					Ystar = self.Ystar
		return Ystar, self.Ystar_info

	def getCorrectedExpression(self,rand_eff_ids=None):
		"""
		Args:
			rand_eff_ids:	index of the random effects that are to consider for the correction
		Returns:
			Ycorr:		corrected expression levels
		"""
		assert self.var is not None, 'scLVM:: use varianceDecomposition method before'

		# check rand_eff_ids
		if rand_eff_ids==None:			rand_eff_ids=range(len(self.Ystar))
		elif type(rand_eff_ids)!=list:	rand_eff_ids=[rand_eff_ids]

		# loop on random effect to consider and correct
		#predicitive means were calculated for standardised expression
		Ystd = self.Y-self.Y.mean(0)
		Ystd/= self.Y.std(0)
		Ycorr = Ystd[:,self.Ystar_info['gene_idx']]#copy.copy(self.Y[:,self.Ystar_info['gene_idx']])
		for i in rand_eff_ids:
			Ycorr -= self.Ystar[i]
		Ycorr*=self.Y[:,self.Ystar_info['gene_idx']].std(0) #bring back to original scale
		Ycorr+=self.Y[:,self.Ystar_info['gene_idx']].mean(0)
		return Ycorr

	def fitLMM(self,K=None,tech_noise=None,idx=None,i0=None,i1=None,verbose=False):
		"""
		Args:
			K:				list of random effects to be considered in the analysis
							if K is none, it does not consider any random effect
			idx:			indices of the genes to be considered in the analysis
			i0:				gene index from which the anlysis starts
			i1:				gene index to which the analysis stops
			verbose:		if True, print progresses
		Returns:
			pv:				matrix of pvalues
			beta:			matrix of correlations
			info:			dictionary annotates pv and beta rows and columns, containing
							gene_idx_row:	index of the genes in rows
							conv:		boolean vetor marking genes for which variance decomposition has converged
							gene_row:   annotate rows of matrices
		"""
		assert self.var is not None, 'scLVM:: when multiple hidden factors are considered, varianceDecomposition decomposition must be used prior to this method'
#		print QTL

		if idx==None:
			if i0==None or i1==None:
				i0 = 0; i1 = self.G
			idx = SP.arange(i0,i1)
		elif type(idx)!=SP.ndarray:
			idx = SP.array([idx])

		if K is not None and type(K)!=list:	K = [K]

		lmm_params = {'covs':SP.ones([self.N,1]),'NumIntervalsDeltaAlt':100,'NumIntervalsDelta0':100,'searchDelta':True}

		Ystd = self.Y-self.Y.mean(0)
		Ystd/= self.Y.std(0)

		beta   = SP.zeros((idx.shape[0],self.G))
		pv	 = SP.zeros((idx.shape[0],self.G))
		geneID = SP.zeros(idx.shape[0],dtype=str)
		count  = 0
		var = self.var/self.var.sum(1)[:,SP.newaxis] 
		for ids in idx:
			if verbose:
				print '.. fitting gene %d'%ids
			# extract a single gene
			if K is not None:
				if len(K)>1:
					if self.var_info['conv'][count]==True:
						_K = SP.sum([var[count,i]*K[i] for i in range(len(K))],0)
						_K/= _K.diagonal().mean()
					else:
						_K = None
				else:
					_K = K[0]
			else:
				_K = None
			lm = QTL.test_lmm(Ystd,Ystd[:,ids:ids+1],K=_K,verbose=False,**lmm_params)
			pv[count,:]   = lm.getPv()[0,:]
			beta[count,:] = lm.getBetaSNP()[0,:]
			if self.geneID is not None:   geneID[count] = self.geneID[ids]
			count+=1

		info = {'conv':self.var_info['conv'],'gene_idx_row':idx}
		if geneID is not None:	info['gene_row'] = geneID

		return pv, beta, info
		
		
			

