""" Generation of heatmaps plotting SPC index for N = 3000
Author 
------
Roman Baravalle <romanbaravalle@gmail.com>

Version
-------
0.1

Date
----
12-7-2021


"""


from scipy.ndimage.filters import gaussian_filter
import pandas

from ImportingPackages import *
from SimulationParametersCritical import *
from AuxiliarFunctions import *
from NeuronModelsVarDrive import *
from SetupModels import *
#execfile("ImportingPackages.py")
#execfile("SimulationParametersCritical.py")
#execfile("AuxiliarFunctions.py")
#execfile("NeuronModels.py")
#execfile("SetupModels.py")

SavingDirectory = os.getcwd()
WorkingDirectory = os.path.dirname(os.path.normpath(SavingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = WorkingDirectory + CSVFolder

NAux = [3000]
alpha = 1
mode = "cluster"
clusters = 1
nu0 = 17
SimuIdentifier = "Raster_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"
kwidth=25.#20.
ksigma=3.#2.

#binedges = np.linspace(1500,1800,81)#1500,1800,81
binedges = np.linspace(1500,1800,81)#1500,1800,81

## Find and open the files
os.chdir(WorkingDirectory)
ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
	J1 = re.search('J=(.+?)_', file1).group(1)
	sigma1 = re.search('sigma=(.+?)_', file1).group(1)
  	Spiketimes = np.loadtxt(file1,skiprows=1,usecols=0)
	NeuronsSpiking = np.unique(np.loadtxt(file1,skiprows=1,usecols=1))

	if np.shape(Spiketimes)[0]>0:		
		STH, binedges1 = np.histogram(Spiketimes, bins=binedges) #Spike time histogram with 1 ms bin. ToDo: remove empty rows -non-spiking neurons-
		if STH is None: continue
		SPC, R2, mean_Tnet = compute_SPC_R2_meanFR(STH,kwidth=kwidth,ksigma=ksigma)
		if R2 is None: R2 = -1.
		if mean_Tnet is None: mean_Tnet = -1.
		if SPC is None: SPC = -1.
		ResultsList = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float), SPC/NAux[0], R2, mean_Tnet]
	else:
		ResultsList = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float), -1., -1., -1.]
	ListsOfLists.append(ResultsList)

DataMatrix = np.row_stack(ListsOfLists)

NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]

HeatMapSPC = 0.1*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapR2 = 0.1*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapTnet = 0.1*np.ones((NumberOfSigmas,NumberOfJs))

indexJ = 0
for J in JAux:
	indexSigma = 0
	CV = []
	for sigma in sigmaAux:
		IndexJ = np.argwhere(abs(DataMatrix[:,0]-J)<0.01)
		IndexSigma = np.argwhere(abs(DataMatrix[:,1]-sigma)<0.0001)
		Indexes = np.intersect1d(IndexJ,IndexSigma)
		if Indexes.size is not 0:		
			HeatMapSPC[indexSigma,indexJ] = DataMatrix[Indexes,2]
			HeatMapR2[indexSigma,indexJ] = DataMatrix[Indexes,3]
			HeatMapTnet[indexSigma,indexJ] = 1/DataMatrix[Indexes,4]*1000.
		indexSigma += 1
	indexJ += 1

dataframe = pandas.DataFrame(HeatMapSPC,index=sigmaAux,columns=JAux)
mask = dataframe == -100.
dataframe[mask] = np.nan
dataframe.index = dataframe.index.to_series().apply(lambda x: np.round(x,4))

dataframe.to_csv(SavingDirectory+"/SPC{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
