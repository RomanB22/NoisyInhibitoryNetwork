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

print(CSVFolder)

CeilingValue = 8 #Superior value for a doublet. Minimum will be 1 ms, due to the synaptic delay time

NAux = [800]
alpha = 1
mode = "cluster"
clusters = 1
nu0 = 17
SimuIdentifier = "ISI_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

def elements(array):
	return array.ndim and array.size

#binedges = np.linspace(2900,3200,301)#np.linspace(1600,1900,301)#

SimuIdentifier = "ISI_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

## Find and open the files
os.chdir(WorkingDirectory) #Change working directory

ListsOfLists = []
NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
	J1 = re.search('J=(.+?)_', file1).group(1)
	sigma1 = re.search('sigma=(.+?)_', file1).group(1)
  	ISI =  np.loadtxt(file1,skiprows=1)
	Doublets =ISI[ISI<CeilingValue]
	if elements(Doublets)>0:
		Results = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float),len(Doublets),len(ISI),1.0*len(Doublets)/len(ISI), np.mean(Doublets),np.std(Doublets)]
	else:
		Results = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float),0,len(ISI),0, 0, 0]

	ListsOfLists.append(Results)

DataMatrix = np.row_stack(ListsOfLists)

NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]
HeatMapAmountDoublets = np.zeros((NumberOfSigmas,NumberOfJs))
HeatMapAmountISI = np.zeros((NumberOfSigmas,NumberOfJs))
HeatMapProportion = np.zeros((NumberOfSigmas,NumberOfJs))
HeatMapMean = np.zeros((NumberOfSigmas,NumberOfJs))
HeatMapSTD = np.zeros((NumberOfSigmas,NumberOfJs))

indexJ = 0
for J in JAux:
	indexSigma = 0
	for sigma in sigmaAux:
		IndexJ = np.argwhere(abs(DataMatrix[:,0]-J)<0.1)
		IndexSigma = np.argwhere(abs(DataMatrix[:,1]-sigma)<0.0001)
		Indexes = np.intersect1d(IndexJ,IndexSigma)
		if Indexes.size is not 0:
			HeatMapAmountDoublets[indexSigma,indexJ] = DataMatrix[Indexes,2][0]
			HeatMapAmountISI[indexSigma,indexJ] = DataMatrix[Indexes,3][0]
			HeatMapProportion[indexSigma,indexJ] = DataMatrix[Indexes,4][0]
			HeatMapMean[indexSigma,indexJ] = DataMatrix[Indexes,5][0]
			HeatMapSTD[indexSigma,indexJ] = DataMatrix[Indexes,6][0]
		indexSigma += 1
	indexJ += 1
dataframeAD = pandas.DataFrame(HeatMapAmountDoublets,index=sigmaAux,columns=JAux)
dataframeAI = pandas.DataFrame(HeatMapAmountISI,index=sigmaAux,columns=JAux)
dataframePD = pandas.DataFrame(HeatMapProportion,index=sigmaAux,columns=JAux)
dataframeMD = pandas.DataFrame(HeatMapMean,index=sigmaAux,columns=JAux)
dataframeSD = pandas.DataFrame(HeatMapSTD,index=sigmaAux,columns=JAux)


dataframeAD.to_csv(SavingDirectory +"/DoubletsAmount_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1]))
dataframeAI.to_csv(SavingDirectory +"/DoubletsAmountISI_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1]))
dataframePD.to_csv(SavingDirectory +"/DoubletsProportion_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1]))
dataframeMD.to_csv(SavingDirectory +"/DoubletsMean_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1]))
dataframeSD.to_csv(SavingDirectory +"/DoubletsStandardDeviation_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1]))
