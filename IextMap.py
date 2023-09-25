""" Generation of heatmaps plotting Total input current for N = 800
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
nu0 = 30
SimuIdentifier = "Summary_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

## Find and open the files
os.chdir(WorkingDirectory)
ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
	with open(file1,'r') as f:
		FullList = []
		for key,group in it.groupby(f,lambda line: line.startswith('Simulation Summary') or \
		line.startswith(' alpha') or line.startswith(' Connection') or\
		line.startswith(' Avg. Individual') or line.startswith(' Synchrony') or line.startswith(' InitialConditions') or line.startswith(' Avg.') or line.startswith(' Binder')):
			if not key:
				List = []
				for line in group:
					Number = re.findall(r'\d+.\d+', line) # find number of digits through regular expression
					if np.size(Number) == 0:
						Number = [0.]
					List.append(float(Number[0]))
				FullList.append(List)
		ResultsList = list(it.chain.from_iterable(FullList))# ResultsList = [I0, J, sigma, N]
		ListsOfLists.append(ResultsList)
"""
ListsOfLists has the four values of [I0, J, sigma, N] for each one of the simulations
So ListsOfLists[0] shows the results for the first simulation
   ListsOfLists[indexSimulation][0] = I0 for simulation = IndexSimulation
   ListsOfLists[indexSimulation][1] = J for simulation = IndexSimulation
   ListsOfLists[indexSimulation][2] = sigma for simulation = IndexSimulation
   ListsOfLists[indexSimulation][3] = N for simulation = IndexSimulation
"""
DataMatrix = np.row_stack(ListsOfLists)


#############

NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]

HeatMapI0 = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapIrecMean = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapIrecMax = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapIrecMin = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapPermanenceSubSNP = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapPermanenceBist = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapPermanenceSupraHopf = -100.*np.ones((NumberOfSigmas,NumberOfJs))
HeatMapPermanenceSupraSN = -100.*np.ones((NumberOfSigmas,NumberOfJs))

indexJ = 0
for J in JAux:
	indexSigma = 0
	for sigma in sigmaAux:
		IndexJ = np.argwhere(abs(DataMatrix[:,1]-J)<0.1)
		IndexSigma = np.argwhere(abs(DataMatrix[:,2]-sigma)<0.0001)
		Indexes = np.intersect1d(IndexJ,IndexSigma)
		if Indexes.size is not 0:
			I0 = DataMatrix[Indexes,0][0]
			HeatMapI0[indexSigma,indexJ] = I0

		indexSigma += 1
	indexJ += 1

dataframeI0 = pandas.DataFrame(HeatMapI0,index=sigmaAux,columns=JAux)
dataframeI0.index = dataframeI0.index.to_series().apply(lambda x: np.round(x,4))

dataframeI0.to_csv(SavingDirectory+"/ExtCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))

print("Early end"); quit()
