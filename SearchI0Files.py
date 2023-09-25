import numpy as np
import glob, os

sigmaAux = np.sort(np.concatenate((np.around(np.logspace(-2,1,41,endpoint=True),4),np.linspace(0.011,10.011,41, endpoint=True))))

JAux = np.sort(np.concatenate((np.around(np.logspace(0,2,21,endpoint=True),1),np.linspace(3,153,31, endpoint=True))))

Folders = {0:'I0_IzhiTyI/', 1:'I0_IzhiTyII/', 2:'I0_LIF/', 3:'I0_AdExTyI/', 4:'I0_AdExTyIIAH/', 5:'I0_AdExTyIISN/', 6:'I0_Guillem/'}

Folder = Folders[4]
N=800
alpha = 1
mode = "cluster"
clusters = 1
nu = 30

Matrix = Folder + 'MatrixI0JSigma_' + Folder[3:-1]+ str(nu) + 'Hz_N' + str(N)
JNotCalculated = []
sigmaNotCalculated = []
MatrixI0 = np.loadtxt(Matrix);
for J in JAux:
	rows = np.argwhere(abs(MatrixI0[:,1]-J)<1e-3)
	for sigma in sigmaAux:
		column = np.argwhere(abs(MatrixI0[rows,2]-sigma)<1e-4)
		if len(column)==0:
			JNotCalculated.append(J)
			sigmaNotCalculated.append(sigma)
print(np.transpose([[JNotCalculated],[sigmaNotCalculated]]));
print(JNotCalculated,sigmaNotCalculated);
print(len(JNotCalculated))

