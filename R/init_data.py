import sys
import scipy as SP
import pylab as PL
from matplotlib import cm
import h5py
limix_path = '/Users/florian/Code/python_code/limix-0.6.4/build/release.darwin/interfaces/python'
sclvm_path = '/Users/florian/Code/python_code/scPy/scLVM/'
sys.path.append(limix_path)
sys.path.append(sclvm_path)
#import scLVM
sys.path.append('./../scLVM')
from scLVM import scLVM


