""" Generation of heatmaps plotting network frequency for N = 800
Author 
------
Roman Baravalle <romanbaravalle@gmail.com>

Version
-------
0.1

Date
----
5-31-2022
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

#sigmaAux = np.linspace(0.011,10.011,41, endpoint=True)
#JAux = np.concatenate(([1],np.linspace(3,153,31, endpoint=True)))

alpha = 1
mode = "cluster"
clusters = 1
nu0 = 30
NAux = [800]

## Find and open the files
os.chdir(WorkingDirectory)

##############
SimuIdentifier = "AvMembVolt*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
	J1 = re.search('J=(.+?)_', file1).group(1)
	sigma1 = re.search('sigma=(.+?)_', file1).group(1)
	Time = np.loadtxt(file1,skiprows=0,usecols=0)
	fs = 1./(Time[1]-Time[0])*1000	
	VoltMean = np.loadtxt(file1,skiprows=0,usecols=1)
	
	ft = np.fft.rfft(VoltMean)
	freqs = np.fft.rfftfreq(len(VoltMean), (Time[1]-Time[0])/1000) # Get frequency axis from the time axis
	mags = abs(ft) # We don't care about the phase information here

	#plt.plot(freqs[:70],mags[:70])
	#plt.savefig("../SimulationAnalysis/AutoCorrelation/Freq{N}_J={J}_sigma={sigma}_{Folder}.png".format(N=NAux[0],J=J1,sigma=sigma1,Folder=CSVFolder[1:-1]))
	#plt.close()

	inflection = np.diff(np.sign(np.diff(mags)))
	peaks = (inflection < 0).nonzero()[0] + 1
	peak = peaks[mags[peaks].argmax()]
	

	if np.shape(peaks)[0]==0:
		signal_freq = -50.
	else:
		signal_freq = freqs[peak] # Gives 0.05

	ResultsList = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float),signal_freq]
	ListsOfLists.append(ResultsList)
DataMatrix = np.row_stack(ListsOfLists)
#############
NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]

HeatMapCrossCorr = -50.*np.ones((NumberOfSigmas,NumberOfJs))

indexJ = 0
for J in JAux:
	indexSigma = 0
	for sigma in sigmaAux:
		IndexJ = np.argwhere(abs(DataMatrix[:,0]-J)<0.1)
		IndexSigma = np.argwhere(abs(DataMatrix[:,1]-sigma)<0.0001)
		Indexes = np.intersect1d(IndexJ,IndexSigma)

		if Indexes.size is not 0:
			Freq = DataMatrix[Indexes,2][0]
			HeatMapCrossCorr[indexSigma,indexJ] = Freq
		indexSigma += 1
	indexJ += 1
 
dataframe = pandas.DataFrame(HeatMapCrossCorr,index=sigmaAux,columns=JAux)

dataframe.to_csv(SavingDirectory +"/NetwFreqMap"+str(nu0)+"CSV_"+CSVFolder[1:-1]+"_N_"+str(NAux[0])+".csv")
