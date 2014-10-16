import scipy as SP
from gp_clvm import gpCLVM
import pdb
import copy
import pylab as PL
PL.ion()

if __name__ == "__main__":

	N = 100 # number of cells
	G = 500	# number of genes
	k0 = 1  # dimension of the known hidden factor (cell cycle)
	k  = 1  # dimension of the hidden factor to infer (t-cells differentiation)

	# generate latent factors
	X0 = SP.randn(N,k0) # known hidden factor (cell cycle)
	X_true = SP.randn(N,k)  # unknown latelt factor

	var = SP.ones(3)/3.

	# generate phenotype
	B0 = SP.randn(k0,G)
	B  = SP.randn(k,G)
	Y_X0   = SP.dot(X0,B0)
	Y_X0  *= SP.sqrt(var[0]/Y_X0.var(0).mean())
	Y_X	= SP.dot(X_true,B) 
	Y_X   *= SP.sqrt(var[1]/Y_X.var(0).mean())
	Y_nois = SP.randn(N,G) 
	Y_nois*= SP.sqrt(var[2]/Y_nois.var(0).mean())
	Y = Y_X0+Y_X+Y_nois
	Y-= Y.mean(0)
	Y/= Y.std(0)

	pdb.set_trace()

	# initialization method
	init_method = 'null' #['fast','regressOut','random','null']
	Ycc = Y_X0+1e-2*SP.randn(N,G)	# cell cycle contributions for regressOut
	X_init = X_true+1e-2*SP.randn(N,k) # initial Xs for regressOut
	varXX = 0.3 # variance explained by the normalized XX.T
	varX0X0 = 0.3 # variance explained by known factor
	nois = 0.7 # variance explained by noise

	if 1:
		"""
		gp with interaction term
		"""
		gp = gpCLVM(Y=Y,X0=X0,k=k)
		params0 = gp.initParams(method=init_method,Ycc=Ycc,X=X_init,varXX=varXX,varX0X0=varX0X0,nois=nois)
		# gp.fix_a(flag=True)
		conv = gp.optimize(params0)

		X = gp.getX() # inferred hidden factor
		K = gp.getK() # matrix from inferred inferred factor
		Ki = gp.getKi() # interaction matrix
		var = gp.getVarianceComps() # variance components

		if 1:
			""" plot """
			PL.subplot(2,2,1)
			PL.plot(X_true.ravel(),X.ravel(),'.k')
			PL.xlabel('true X')
			PL.ylabel('fitted X')
			plt = PL.subplot(2,2,2)
			varV = SP.array([var[key] for key in var.keys()])
			varV/= varV.sum()
			PL.bar(SP.arange(varV.shape[0]),varV)
			plt.set_xticks(SP.arange(varV.shape[0])+0.5)
			plt.set_xticklabels(var.keys())
			PL.savefig('./figures/X.pdf')
			pdb.set_trace()
			PL.clf()

		pdb.set_trace()

	if 1:
		"""
		gp without interaction term
		"""
		_gp = gpCLVM(Y=Y,X0=X0,k=k,interaction=False)
		_params0 = _gp.initParams(method=init_method,Ycc=Ycc,X=X_init,varXX=varXX)
		conv = _gp.optimize(_params0)

		_X = _gp.getX() # inferred hidden factor
		_K = _gp.getK() # matrix from inferred inferred factor
		_var = _gp.getVarianceComps() # variance components

		if 1:
			""" plot """
			PL.subplot(2,2,1)
			PL.plot(X_true.ravel(),_X.ravel(),'.k')
			PL.xlabel('true X')
			PL.ylabel('fitted X')
			plt = PL.subplot(2,2,2)
			_varV = SP.array([_var[key] for key in _var.keys()])
			_varV/= _varV.sum()
			PL.bar(SP.arange(_varV.shape[0]),_varV)
			plt.set_xticks(SP.arange(_varV.shape[0])+0.5)
			plt.set_xticklabels(_var.keys())
			PL.savefig('./figures/X_noint.pdf')

		pdb.set_trace()


