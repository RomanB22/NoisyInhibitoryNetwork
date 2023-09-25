"""
File to load all the I0 files, and saving all the info as a matrix in a single file

"""
import glob, os
import numpy as np
execfile("FoldersName.py")

ChosenModel = 0 #Models = {0:'IzhiTypeI', 1:'IzhiTypeII', 2:'LIFModel', 3:'AdExModel', 4:'GuillemModel'}
ChosenSubmodelAdEx = 2 #SubModelsAdEx = {0:'TypeI', 1:'TypeIISaddleNode', 2:'TypeIIAndronovHopf'}
ThetaInput = False

if ChosenModel==0:
	I0Folder = I0FolderIzhiTypeI
elif ChosenModel==1:
	I0Folder = I0FolderIzhiTypeII
elif ChosenModel==2:
	I0Folder = I0FolderLIF
elif ChosenModel==3:
	if ChosenSubmodelAdEx==0:
		I0Folder = I0FolderAdExTypeI
	elif ChosenSubmodelAdEx==1:
		I0Folder = I0FolderAdExTypeIISN
	elif ChosenSubmodelAdEx==2:
		I0Folder = I0FolderAdExTypeIIAndrHopf
elif ChosenModel==4:
	I0Folder = I0FolderGuillem
elif ChosenModel==5:
	I0Folder = I0FolderIzhiTypeIIMod

WorkingDirectory = os.getcwd() + I0Folder

N=3000
alpha = 1
mode = "cluster"
clusters = 1
nu = 30

SimuIdentifier = '*_N=' + str(N) + '*_alpha=' + str(alpha) + '*_nu0=' + str(nu) + '*_InitCond=' + str(mode) + str(clusters) + '*.csv'

print(SimuIdentifier,WorkingDirectory)

## Find and open the files
os.chdir(WorkingDirectory)
ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))
filename = "MatrixI0JSigma" + I0Folder[3:-1] + str(nu) + "Hz" + "_N" + str(N)
if os.path.exists(filename):
	MatrixDataNP = np.loadtxt(filename)
else:
	MatrixDataNP=[];
	print("No Matrix");
MatrixData = []

for file1 in glob.glob(SimuIdentifier):
	A = np.loadtxt(file1)
	MatrixData.append(A)
if len(MatrixData)!=0:
	MatrixData = np.row_stack(MatrixData)
	if len(MatrixDataNP)!=0:	
		MatrixData = np.concatenate((MatrixDataNP,MatrixData))	
else:
	MatrixData = MatrixDataNP
MatrixData1, Index = np.unique(MatrixData[:,1:3], axis=0, return_index = True)

print(np.shape(MatrixData[Index,:]))
np.savetxt(filename,MatrixData[Index,:],fmt="%2.6f")
