import matplotlib.pyplot as plt
import numpy as np
import os, glob
import re
from AuxiliarFunctions import compute_SPC_R2_meanFR

#sigmaAux = np.around(np.linspace(0.011,10.011,51, endpoint=True),4)
#JAux = [1.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.]
#I0Aux = [0.,5.,10.,15.,20.,25.,30.,35.,40.]
#JAux = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
#I0Aux= [10.,20.,30.,40.,50.,60.,70.]#[0.,10.,20.,30.,40.,50.,60.,70.]

#Conductance based
sigmaAux = np.concatenate ( (np.around(np.linspace(0.011,10.011,51, endpoint=True),4) , np.around(np.linspace(10.011,15.011,11, endpoint=True),4)) )#np.around(np.linspace(0.011,10.011,51, endpoint=True),4) 
JAux = [1.,20.,40.,60.,80.,100.,120.,140.,160.,180.,250.,300.]
I0Aux= [60.,75.]#[30.,40.,45.,50.,55.,60.,65.,70.,75.,80.]

#JAux = [1.6,79.4]
#I0Aux= [0.,5.,10.,16.,20.,25.,30.,2.03125,3.59375,3.984375,4.2773,7.685547,15.29297]


WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_GuillemVarDriveHyper/"
#WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_IzhiTyIIVarDrive/"

NAux = [3000]
alpha = 1
mode = "cluster"
clusters = 1

SimuIdentifier = "Raster_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*I0=75.*_InitCond=' + str(mode) + str(clusters) + "*-75.0.csv"

SimuIdentifierFullMon = "FullMonitor_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*I0=75.*_InitCond=' + str(mode) + str(clusters) + "*-75.0.csv"  

kwidth=25.#20.#20.25.15.
ksigma=3.#2.#2.3.1.5
binedges = np.linspace(1500,1900,71)#np.linspace(1100,1500,91)

## Find and open the files
os.chdir(WorkingDirectory)

ListsOfLists = []

for file1 in glob.glob(SimuIdentifierFullMon):
    print(file1)
    J1 = float(re.search('J=(.+?)_', file1).group(1))
    sigma1 = float(re.search('sigma=(.+?)_', file1).group(1))
    I01 = float(re.search('I0=(.+?)_', file1).group(1))

    FullCurrent = np.loadtxt(file1,skiprows=1,usecols=2,ndmin=1)
    ResultsList = [J1, I01,  sigma1, np.max(FullCurrent), np.mean(FullCurrent)]
    ListsOfLists.append(ResultsList)
DataMatrix = np.row_stack(ListsOfLists)

np.savetxt("/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes/DataMatrixCurrentGuillemHyper_I075.txt",DataMatrix,fmt="%3.8f")



ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
    J1 = float(re.search('J=(.+?)_', file1).group(1))
    sigma1 = float(re.search('sigma=(.+?)_', file1).group(1))	
    I01 = float(re.search('I0=(.+?)_', file1).group(1))

    Spiketimes = np.loadtxt(file1,skiprows=1,usecols=0,ndmin=1)
    NeuronsSpiking = np.unique(np.loadtxt(file1,skiprows=1,usecols=1,ndmin=1))

    if np.shape(Spiketimes)[0]>0:
        STH, binedges1 = np.histogram(Spiketimes, bins=binedges) #Spike time histogram with 1 ms bin. ToDo: remove empty rows -non-spiking neurons-
	
        if STH is None: continue
        SPC, R2, mean_Tnet = compute_SPC_R2_meanFR(STH,kwidth=kwidth,ksigma=ksigma)
        if R2 is None: R2 = -1.
        if mean_Tnet is None: mean_Tnet = -1.
        if SPC is None: SPC = -1.
        ResultsList = [J1, I01,  sigma1, SPC/NAux[0], R2, mean_Tnet,len(NeuronsSpiking)]
    else:
        ResultsList = [J1, I01,  sigma1, -1., -1., -1.,len(NeuronsSpiking)]
    ListsOfLists.append(ResultsList)

DataMatrix = np.row_stack(ListsOfLists)
np.savetxt("/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes/DataMatrixSPCGuillemHyper_I075.txt",DataMatrix,fmt="%3.8f")

#J1, I01,  sigma1, SPC, R2, mean_Tnet = np.loadtxt("DataMatrixSPC.txt",skiprows=1,unpack=True)
