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

"""Conditions CLVM code.
This is currently experimental and not supported
"""

import pdb
import scipy as SP
import scipy.linalg as LA
import copy
import pdb

import sys
sys.path.append('./..')
from utils.misc import regressOut

# import limix

import limix

class gpCLVM(object):
	"""
	Class for conditional gplvm
	"""

	def __init__(self,Y=None,X0=None,k=1,standardize=False,interaction=True):
		"""
		Y:		data [NxG]
		X0:		known latent factors [Nxk0]
		k:		number of latent factors to infer
		"""

		assert Y!=None, 'gpCLVM:: set Y!'
		assert X0!=None, 'gpCLVM:: set X0!'

		self.cache = {}
		self.interaction = interaction

		# read data
		self._setY(Y,standardize=standardize)
		self._setX0(X0)
		self.k = k

		# covariance for known latex factor
		self.C0  = limix.CFixedCF(self.K0)

		# covariance for unknow latent factor
		self.C  = limix.CProductCF()
		self.Ca = limix.CLowRankCF(self.N,self.k)
		self.C.addCovariance(self.Ca)
		if self.interaction==True:
			self.Cb1 = limix.CFixedCF(SP.ones((self.N,self.N)))
			self.Cb1.setParamMask(SP.zeros(1))
			self.Cb2 = limix.CFixedCF(self.K0)
			self.Cb  = limix.CSumCF()
			self.Cb.addCovariance(self.Cb1)
			self.Cb.addCovariance(self.Cb2)
			self.C.addCovariance(self.Cb)

		# total covariance
		covar  = limix.CSumCF()
		covar.addCovariance(self.C0)
		covar.addCovariance(self.C)

		# likelihood
		self.ll = limix.CLikNormalIso()

		# init GP and hyper params
		self.gp = limix.CGPbase(covar,self.ll)
		self.gp.setY(self.Y)

	def _setY(self,Y,standardize=False):
		""" set phenotype """
		Y -= Y.mean(0)
		if standardize:
			Y /= Y.std(0)
		self.Y  = Y
		self.N,self.G=Y.shape

	def _setX0(self,X0):
		""" set X0 """
		assert X0.shape[0]==self.N, 'gpCLVM:: dimension dismatch'
		self.X0 = X0
		self.k0 = X0.shape[0]
		self.K0  = SP.dot(self.X0,self.X0.T)
		self.K0 /= self.K0.diagonal().mean()

	def initParams(self,method='fast',Ycc=None,X=None,varXX=None,varX0X0=None,nois=None):
		"""
		This takes care of parameter initialization
		"""
		if method=='fast':
			rv = self._initParams_fast()
		elif method=='regressOut':
			assert Ycc!=None, 'provide Ycc'
			assert X!=None, 'provide X'
			assert varXX!=None, 'provide varXX'
			rv = self._initParams_regressOut(Ycc,X,varXX)
		elif method=='random':
			rv = self._initParams_random()
		elif method=='null':
			assert varX0X0!=None, 'provide varX0X0'
			assert nois!=None, 'provide nois'
			rv = self._initParams_null(varX0X0,nois)
		return rv

	def optimize(self,params0):
		""" initialize """
		self.gp.setParams(params0)
		self.cache['params0'] = params0 
		self.cache['lml0'] = self.gp.LML()
		self.cache['lmlGrad0'] = self.gp.LMLgrad()
		""" optimize """
		gpopt = limix.CGPopt(self.gp)
		conv = gpopt.opt()
		""" store stuff """
		self.cache['lml'] = self.gp.LML()
		self.cache['lmlGrad'] = self.gp.LMLgrad()
		self.cache['params'] = self.gp.getParams()
		self.cache['X'] = self.Ca.getParams().reshape((self.N,self.k),order='F')
		self.cache['K'] = self.Ca.K()
		if self.interaction:
			self.cache['Ki'] = self.Ca.K()*self.Cb2.K()
		else:
			self.cache['Ki'] = None
		self.cache['var'] = {}
		self.cache['var']['K0'] = self.C0.getParams()[0]**2
		self.cache['var']['K']  = self.cache['K'].diagonal().mean()
		if self.interaction:
			self.cache['var']['Ki'] = self.cache['Ki'].diagonal().mean()
		self.cache['var']['noise'] = self.ll.getParams()[0]**2
		self.cache['K']  /= self.cache['var']['K']
		if self.interaction:
			self.cache['Ki'] /= self.cache['var']['Ki']
		return conv

	def getX(self):
		"""
		return X
		"""
		return self.cache['X']

	def getK(self):
		"""
		return K
		"""
		return self.cache['K']

	def getKi(self):
		"""
		return Ki
		"""
		return self.cache['Ki']
		
	def getVarianceComps(self):
		"""
		return variance compoennts
		"""
		return self.cache['var']


	def _initParams_fast(self):
		""" 
		initialize the gp parameters
			1) project Y on the known factor X0 -> Y0
				average variance of Y0 is used to initialize the variance explained by X0
			2) considers the residual Y1 = Y-Y0 (this equivals to regress out X0)
			3) perform PCA on cov(Y1) and considers the first k PC for initializing X
			4) the variance of all other PCs is used to initialize the noise
			5) the variance explained by interaction is set to a small random number 
		"""
		Xd = LA.pinv(self.X0)
		Y0 = self.X0.dot(Xd.dot(self.Y))
		Y1 = self.Y-Y0
		YY = SP.cov(Y1)
		S,U = LA.eigh(YY)
		X = U[:,-self.k:]*SP.sqrt(S[-self.k:])
		a = SP.array([SP.sqrt(Y0.var(0).mean())])
		b = 1e-3*SP.randn(1)
		c = SP.array([SP.sqrt((YY-SP.dot(X,X.T)).diagonal().mean())])
		# gp hyper params
		params = limix.CGPHyperParams()
		if self.interaction:
			params['covar'] = SP.concatenate([a,X.reshape(self.N*self.k,order='F'),SP.ones(1),b])
		else:
			params['covar'] = SP.concatenate([a,X.reshape(self.N*self.k,order='F')])
		params['lik'] = c
		return params
		

	def _initParams_regressOut(self,Ycc,X,varXX):
		""" 
		initialize the gp parameters
			1) the variance of Kcc as Ycc.var(0).mean()
			2) X with the provided 
			3) variance of interaction (if label is True) will be set to ~0
			4) residual to residual
		"""
		X *= SP.sqrt(varXX/(X**2).mean())
		Y1 = self.Y-Ycc
		a = SP.array([SP.sqrt(Ycc.var(0).mean())])
		b = 1e-3*SP.ones(1)
		c = Y1.var(0).mean()-varXX
		c = SP.maximum(1e-1,c)
		c = SP.array([SP.sqrt(c)])
		# gp hyper params
		params = limix.CGPHyperParams()
		if self.interaction:
			params['covar'] = SP.concatenate([a,X.reshape(self.N*self.k,order='F'),SP.ones(1),b])
		else:
			params['covar'] = SP.concatenate([a,X.reshape(self.N*self.k,order='F')])
		params['lik'] = c
		return params
		

	def _initParams_random(self):
		""" 
		initialize the gp parameters randomly
		"""
		# gp hyper params
		params = limix.CGPHyperParams()
		if self.interaction:	params['covar'] = SP.concatenate([SP.randn(self.N*self.k+1),SP.ones(1),SP.randn(1)])
		else:					params['covar'] = SP.randn(self.N*self.k+1)
		params['lik'] = SP.randn(1)
		return params
		
	def _initParams_null(self,varX0X0,nois):
		""" 
		initialize from null model
		"""
		X  = 1e-3*SP.randn(self.N,self.k)
		a = SP.array([SP.sqrt(varX0X0)])
		b = 1e-3*SP.ones(1)
		c = SP.array([SP.sqrt(nois)])
		# gp hyper params
		params = limix.CGPHyperParams()
		if self.interaction:
			params['covar'] = SP.concatenate([a,X.reshape(self.N*self.k,order='F'),SP.ones(1),b])
		else:
			params['covar'] = SP.concatenate([a,X.reshape(self.N*self.k,order='F')])
		params['lik'] = c
		return params

	def fix_a(self,flag=True):
		"""
		if flag==True:	fix a
		else:			set a free
		"""
		if flag:	self.C0.setParamMask(SP.zeros(1))
		else:		self.C0.setParamMask(SP.ones(1))
		
