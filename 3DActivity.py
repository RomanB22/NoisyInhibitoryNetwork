import matplotlib.pyplot as plt
import numpy as np
import os, glob
import re
from AuxiliarFunctions import compute_SPC_R2_meanFR

sigmaAux = np.unique( np.concatenate ( (np.around(np.linspace(0.011,10.011,51, endpoint=True),4) , np.around(np.linspace(10.011,15.011,11, endpoint=True),4)) ) )

JAux = [1.6,79.4,100.]
I0Aux = np.sort([0.,5.,10.,16.,20.,25.,30.,2.03125,3.59375,3.984375,4.2773,7.685547,15.29297])
#I0Aux = np.sort([5.,10.,16.,20.,25.,30.,2.03125,3.59375,3.984375,4.2773,7.685547,15.29297])

WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_GuillemVarDriveShunt/"

NAux = [3000]
alpha = 1
mode = "cluster"
clusters = 1

SimuIdentifier = "Raster_*J=100.*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"
SimuIdentifierFullMon = "FullMonitor_*J=100.*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"  

kwidth=20
ksigma=2
#binedges = np.linspace(900,1300,101)
binedges = np.linspace(900,1300,151)
## Find and open the files
os.chdir(WorkingDirectory)

ListsOfLists = []

for file1 in glob.glob(SimuIdentifierFullMon):
    J1 = float(re.search('J=(.+?)_', file1).group(1))
    sigma1 = float(re.search('sigma=(.+?)_', file1).group(1))
    I01 = float(re.search('I0=(.+?)_', file1).group(1))

    if np.any(abs( I0Aux-I01 )<1e-1) and np.any(abs( sigmaAux-sigma1 )<1e-3):
        FullCurrent = np.loadtxt(file1,skiprows=1,usecols=3,ndmin=1)
        ResultsList = [I01,  sigma1, np.max(FullCurrent)+10*I01, np.mean(FullCurrent)+10*I01]

        ListsOfLists.append(ResultsList)

DataMatrix = np.row_stack(ListsOfLists)

np.savetxt("/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes/DataMatrixCurrentIzhiTyII_100.txt",DataMatrix,fmt="%1.4f")



ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
    J1 = float(re.search('J=(.+?)_', file1).group(1))
    sigma1 = float(re.search('sigma=(.+?)_', file1).group(1))	
    I01 = float(re.search('I0=(.+?)_', file1).group(1))

    if np.any(abs( I0Aux-I01 )<1e-1) and np.any(abs( sigmaAux-sigma1 )<1e-3):
        Spiketimes = np.loadtxt(file1,skiprows=1,usecols=0,ndmin=1)
        NeuronsSpiking = np.unique(np.loadtxt(file1,skiprows=1,usecols=1,ndmin=1))
        if np.shape(Spiketimes)[0]>0:
            STH, binedges1 = np.histogram(Spiketimes, bins=binedges) #Spike time histogram with 1 ms bin. ToDo: remove empty rows -non-spiking neurons-
	
            if STH is None: continue
            SPC, R2, mean_Tnet = compute_SPC_R2_meanFR(STH,kwidth=kwidth,ksigma=ksigma)
            if R2 is None: R2 = -1.
            if mean_Tnet is None: mean_Tnet = -1.
            if SPC is None: SPC = -1.
            ResultsList = [I01,  sigma1, SPC/NAux[0], R2]
        else:
            ResultsList = [I01,  sigma1, -1., -1.]
        ListsOfLists.append(ResultsList)

DataMatrix = np.row_stack(ListsOfLists)
np.savetxt("/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes/DataMatrixSPCIzhiTyII_100.txt",DataMatrix,fmt="%1.5f")

#J1, I01,  sigma1, SPC, R2, mean_Tnet = np.loadtxt("DataMatrixSPC.txt",skiprows=1,unpack=True)
