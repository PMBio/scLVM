# Copyright(c) 2014, The scLVM developers (Forian Buettner, Paolo Francesco Casale, Oliver Stegle)
#
#Licensed under the Apache License, Version 2.0 (the "License");
#you may not use this file except in compliance with the License.
#You may obtain a copy of the License at
#
#	http://www.apache.org/licenses/LICENSE-2.0
#
#Unless required by applicable law or agreed to in writing, software
#distributed under the License is distributed on an "AS IS" BASIS,
#WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#See the License for the specific language governing permissions and
#limitations under the License.

import sys
sys.path.append('./..')
import limix_legacy

import limix_legacy.deprecated.modules.panama as PANAMA
import limix_legacy.deprecated.modules.varianceDecomposition as VAR
import limix_legacy.deprecated.modules.qtl as QTL    
	
import scipy as SP
import scipy.linalg
import scipy.stats
import pdb
import h5py
import time
import copy
import warnings
import os
from utils.misc import dumpDictHdf5
from utils.misc import PCA 
from utils.misc import warning_on_one_line 
from gp_clvm import gpCLVM

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
		if tech_noise!=None:
			self.set_tech_noise(tech_noise)

	
	def fitFactor(self,idx=None,X0=None,k=1,standardize=False, use_ard=False, interaction=True, initMethod='fast'):
		"""
		Args:
			idx:			index of the genes involved
							(e.g., for capturing cell cycle, index of cell cycle genes)
			X0:			known factor(s) on which to condition on when fitting the GPLVM
			k:				number of latent factors
			standardize:	if True, rescale gene expression by std prior to fitting GPLVM
							(data are always mean-centered)
			use_ard:		use automatic relevance detection (switch off unimportant factors)
		Returns:
			X:				hidden variable
			Kconf:			similarity matrix based on the confounding effect (XX.T)
			varGPLVM:		variance contributions of latent factors and residual biological noise
		"""
		assert idx!=None, 'scLVM:: specify idx'
		if use_ard==True and  k<2:
			warnings.formatwarning = warning_on_one_line
			warnings.warn('when using ARD consider choosing k>1')
		if X0!=None:
			assert use_ard == False, 'scLVM:: when fitting conditional GPLVM, use_ard has to be False'		

		Yconf = self.Y[:,idx]
		Yconf-= Yconf.mean(0)
		Kint = None
		if X0==None:
			#use PANAMA
			# prepare data
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
		else:
			# use gpCLVM
			gp = gpCLVM(Y=Yconf,X0=X0,k=k,standardize=standardize,interaction=interaction)
			params0 = gp.initParams(method=initMethod)
			conv = gp.optimize(params0)
			X = gp.getX()
			Kconf = gp.getK()
			if interaction==True:
				Kint = gp.getKi()
			varGPLVM = gp.getVarianceComps()
		return X,Kconf,Kint,varGPLVM

	def set_tech_noise(self,tech_noise):
		"""
		Args:
			tech_noise:		G vector of technical noise
		"""
		assert tech_noise.shape[0]==self.G, 'scLVM:: tech_noise dimension dismatch'
		
		self.tech_noise = tech_noise
		
	def varianceDecomposition(self,K=None,tech_noise=None,idx=None,i0=None,i1=None,max_iter=10,verbose=False, cache=True):
		"""
		Args:
			K:				list of random effects to be considered in the analysis
			idx:			indices of the genes to be considered in the analysis
			i0:				gene index from which the anlysis starts
			i1:				gene index to which the analysis stops
			max_iter:		maximum number of random restarts
			verbose:		if True, print progresses
		"""

		if tech_noise!=None:		self.set_tech_noise(tech_noise)
		assert self.tech_noise!=None, 'scLVM:: specify technical noise'
		assert K!=None, 'scLVM:: specify K'

		if type(K)!=list:	K = [K]
		for k in K:
			assert k.shape[0]==self.N, 'scLVM:: K dimension dismatch'
			assert k.shape[1]==self.N, 'scLVM:: K dimension dismatch'

		if idx==None:
			if i0==None or i1==None:
				i0 = 0; i1 = self.G
			idx = SP.arange(i0,i1)
		elif type(idx)!=SP.ndarray:
			idx = SP.array(idx)
		idx = SP.intersect1d(SP.array(idx),SP.where(self.Y.std(0)>0)[0]) #only makes sense if gene is expressed in at least one cell
		_G	 = len(idx)
		var	= SP.zeros((_G,len(K)+2))
		_idx   = SP.zeros(_G)
		geneID = SP.zeros(_G,dtype=str)
		conv   = SP.zeros(_G)==1
		Ystar  = [SP.zeros((self.N,_G)) for i in range(len(K))]
		count  = 0
		Yidx = self.Y[:,idx]
		Ystd = Yidx-Yidx.mean(0) 
		Ystd/= Yidx.std(0) #delta optimization might be more efficient
		tech_noise = self.tech_noise[idx]/SP.array(Yidx.std(0))**2

		for ids in range(_G):
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
				continue
			_var = vc.getVarianceComps()[0,:]
			KiY = vc.gp.agetKEffInvYCache().ravel()
			for ki in range(len(K)):
				Ystar[ki][:,count]=_var[ki]*SP.dot(K[ki],KiY)
			var[count,:] = _var
			count+=1;
		if self.geneID!=None:	geneID = SP.array(self.geneID)[idx]
		col_header = ['hidden_%d'%i for i in range(len(K))]
		col_header.append('biol_noise')
		col_header.append('tech_noise')
		col_header = SP.array(col_header)

		#annotate column and rows of var and Ystar
		var_info = {'gene_idx':idx,'col_header':col_header,'conv':conv}
		if geneID!=None:	var_info['geneID'] = SP.array(geneID)
		Ystar_info = {'gene_idx':idx,'conv':conv}
		if geneID!=None:	Ystar_info['geneID'] = SP.array(geneID)

		# cache stuff
		if cache == True:
			self.var   = var
			self.Ystar = Ystar
			self.var_info   = var_info
			self.Ystar_info = Ystar_info
		else:
			return var, var_info

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
		assert self.var!=None, 'scLVM:: use varianceDecomposition method before'
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
		assert self.var!=None, 'scLVM:: use varianceDecomposition method before'
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
		assert self.var!=None, 'scLVM:: use varianceDecomposition method before'

		# check rand_eff_ids
		if rand_eff_ids==None:			rand_eff_ids=range(len(self.Ystar))
		elif type(rand_eff_ids)!=list:	rand_eff_ids=[rand_eff_ids]

		# loop on random effect to consider and correct
		#predicitive means were calculated for standardised expression
		idx = self.Ystar_info['gene_idx']
		Yidx = self.Y[:,idx]
		Ystd = Yidx-Yidx.mean(0)
		Ystd/= Yidx.std(0) #delta optimization might be more efficient
		Ycorr = Ystd#copy.copy(self.Y[:,self.Ystar_info['gene_idx']])
		for i in rand_eff_ids:
			Ycorr -= self.Ystar[i]
		Ycorr*=self.Y[:,self.Ystar_info['gene_idx']].std(0) #bring back to original scale
		Ycorr+=self.Y[:,self.Ystar_info['gene_idx']].mean(0)
		return Ycorr

	def fitLMM(self,expr = None,K=None,tech_noise=None,idx=None,i0=None,i1=None,verbose=False, recalc=True, standardize=True):
		"""
		Args:
			K:				list of random effects to be considered in the analysis
							if K is none, it does not consider any random effect
			expr:				correlations are calculated between the gene expression data (self.Y) and these measures provided in expr. If None, self.Y i sused 	
			idx:
			indices of the genes to be considered in the analysis
			i0:				gene index from which the anlysis starts
			i1:				gene index to which the analysis stops
			verbose:		if True, print progress
			recalc:			if True, re-do variance decomposition
			standardize:		if True, standardize also expression 
		Returns:
			pv:				matrix of pvalues
			beta:			matrix of correlations
			info:			dictionary annotates pv and beta rows and columns, containing
							gene_idx_row:	index of the genes in rows
							conv:		boolean vetor marking genes for which variance decomposition has converged
							gene_row:   annotate rows of matrices
		"""
		

		if idx==None:
			if i0==None or i1==None:
				i0 = 0; i1 = self.G
			idx = SP.arange(i0,i1)
		elif type(idx)!=SP.ndarray:
			idx = SP.array(idx)
		idx = SP.intersect1d(idx,SP.where(self.Y.std(0)>0)[0]) #only makes sense if gene is expressed in at least one cell

		
		if K!=None:
			if type(K)!=list:	K = [K]
			if (recalc==True and len(K)>1) or (recalc==True and self.var==None):
				print 'performing variance decomposition first...'
				var_raw,var_info = self.varianceDecomposition(K=K,idx=idx, cache=False) 
				var = var_raw/var_raw.sum(1)[:,SP.newaxis]
			elif recalc==False and len(K)>1:
				assert self.var!=None, 'scLVM:: when multiple hidden factors are considered, varianceDecomposition decomposition must be used prior to this method'
				warnings.warn('scLVM:: recalc should only be set to False by advanced users: scLVM then assumes that the random effects are the same as those for which the variance decompostion was performed earlier.')
				var_raw = self.var
 				var_info = self.var_info
				var = var_raw/var_raw.sum(1)[:,SP.newaxis]
		
		lmm_params = {'covs':SP.ones([self.N,1]),'NumIntervalsDeltaAlt':100,'NumIntervalsDelta0':100,'searchDelta':True}
				
				
		Yidx = self.Y[:,idx]
		Ystd = Yidx-Yidx.mean(0)
		Ystd/= Yidx.std(0) #delta optimization might be more efficient
		
		if expr==None:
			expr = Ystd		
		elif standardize==True:
			exprStd = expr
			exprStd = expr-expr.mean(0)
			exprStd/= expr.std(0)
			expr = exprStd

		_G1	  = idx.shape[0]
		_G2	 = expr.shape[1]

		geneID = SP.zeros(_G1,dtype=str)
		
		beta   = SP.zeros((_G1,_G2))
		pv	 = SP.zeros((_G1,_G2))
		count  = 0
		
		for ids in range(_G1):
			if verbose:
				print '.. fitting gene %d'%ids
			# extract a single gene
			if K!=None:
				if len(K)>1:
					if var_info['conv'][count]==True:
						_K = SP.sum([var[count,i]*K[i] for i in range(len(K))],0)
						_K/= _K.diagonal().mean()
					else:
						_K = None
				else:
					_K = K[0]
			else:
				_K = None
			lm = QTL.test_lmm(expr,Ystd[:,ids:ids+1],K=_K,verbose=False,**lmm_params)   
			pv[count,:]   = lm.getPv()[0,:]
			beta[count,:] = lm.getBetaSNP()[0,:]
			count+=1

		if self.geneID!=None:   geneID = SP.array(self.geneID)[idx]
		if recalc==True and K!=None  and len(K)>1:	
			info = {'conv':var_info['conv'],'gene_idx_row':idx}
		else:	
			info = {'gene_idx_row':idx}
		if geneID!=None:	info['gene_row'] = geneID

		return pv, beta, info
			
