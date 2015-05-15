import sys
import scipy as SP
import pylab as PL
from matplotlib import cm
import h5py
import warnings
#make sure your paths point to limix and scLVM directories
#limix_path = '/Users/flo/software/limix-master/build/release.darwin/interfaces/python'
#sys.path.append(limix_path)
sys.path.append(sclvm_path)
#from scLVM import scLVM

import limix
from pysrc.core import scLVM
from pysrc.gp_clvm import gpCLVM
#switch off FutureWarnings
warnings.simplefilter(action="ignore", category=FutureWarning) 

