# Helper functions

import sys
import scipy as SP
import h5py
import pdb
import scipy.linalg as LA

def smartAppend(table,name,value):
	"""
	helper function for apppending in a dictionary	
	"""	
	if name not in table.keys():
		table[name] = []
	table[name].append(value)

def dumpDictHdf5(RV,o):
	""" Dump a dictionary where each page is a list or an array """
	for key in RV.keys():
		o.create_dataset(name=key,data=SP.array(RV[key]),chunks=True,compression='gzip')

def smartDumpDictHdf5(RV,o):
	""" Dump a dictionary where each page is a list or an array or still a dictionary (in this case, it iterates)"""
	for key in RV.keys():
		if type(RV[key])==dict:
			g = o.create_group(key)
			smartDumpDictHdf5(RV[key],g)
		else:
			o.create_dataset(name=key,data=SP.array(RV[key]),chunks=True,compression='gzip')

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
    return ' %s:%s: %s:%s' % (filename, lineno, category.__name__, message)

def regressOut(Y,X):
	"""
	regresses out X from Y
	"""
	Xd = LA.pinv(X)
	Y_out = Y-X.dot(Xd.dot(Y))
	return Y_out

def PCA(Y, components):
	"""run PCA, retrieving the first (components) principle components
	return [s0, eig, w0]
	s0: factors
	w0: weights
	"""
	sv = LA.svd(Y, full_matrices=0);
	[s0, w0] = [sv[0][:, 0:components], SP.dot(SP.diag(sv[1]), sv[2]).T[:, 0:components]]
	v = s0.std(axis=0)
	s0 /= v;
	w0 *= v;
	return [s0, w0]
