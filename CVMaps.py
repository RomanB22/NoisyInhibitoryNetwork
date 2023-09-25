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
from scipy.stats import variation 

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

NAux = [800]
alpha = 1
mode = "random"
clusters = 1
nu0 = 17
print(SavingDirectory,nu0)
SimuIdentifier = "ISI_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

binedges = np.linspace(1500,1800,151)

## Find and open the files
os.chdir(WorkingDirectory)
ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
	J1 = re.search('J=(.+?)_', file1).group(1)
	sigma1 = re.search('sigma=(.+?)_', file1).group(1)
  	ISI = np.loadtxt(file1,skiprows=1,usecols=0)
	if np.shape(ISI)[0]>0:	
		CV = variation(ISI)	
		ResultsList = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float), CV]
	else:
		ResultsList = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float), -1.]
	ListsOfLists.append(ResultsList)

DataMatrix = np.row_stack(ListsOfLists)

NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]

HeatMapSPC = -0.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapR2 = -0.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapTnet = -0.*np.ones((NumberOfSigmas,NumberOfJs))

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
		indexSigma += 1
	indexJ += 1

dataframe = pandas.DataFrame(HeatMapSPC,index=sigmaAux,columns=JAux)
mask = dataframe == -100.
dataframe[mask] = np.nan
dataframe.index = dataframe.index.to_series().apply(lambda x: np.round(x,4))

dataframe.to_csv(SavingDirectory+"/CVRandom{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
