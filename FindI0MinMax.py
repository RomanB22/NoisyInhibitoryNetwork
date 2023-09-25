import numpy as np
import glob, os, sys

Folders = {0:'I0_IzhiTyI/', 1:'I0_IzhiTyII/', 2:'I0_LIF/', 3:'I0_AdExTyI/', 4:'I0_AdExTyIIAH/', 5:'I0_AdExTyIISN/', 6:'I0_Guillem/'}

Folder = Folders[6]
N=sys.argv[1]
alpha = 1
mode = "cluster"
clusters = 1
nu = sys.argv[2]

Matrix = Folder + 'MatrixI0JSigma_' + Folder[3:-1]+ str(nu) + 'Hz_N' + str(N)
MatrixI0 = np.loadtxt(Matrix);
print(np.min(MatrixI0[:,0]),np.max(MatrixI0[:,0]))
