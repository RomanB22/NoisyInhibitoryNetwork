import matplotlib.pyplot as plt
import numpy as np
import os, glob
import re
from scipy import signal

#sigmaAux = np.around(np.linspace(0.011,10.011,51, endpoint=True),4)
#JAux = [1.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.]
#I0Aux = [0.,5.,10.,15.,20.,25.,30.,35.,40.]
#JAux = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
#I0Aux= [10.,20.,30.,40.,50.,60.,70.]#[0.,10.,20.,30.,40.,50.,60.,70.]
sigmaAux = np.around(np.linspace(0.011,30.011,51, endpoint=True),4)
JAux = [0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
I0Aux= [10.,20.,30.,40.,50.,60.,70.,80]

WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_GuillemVarDriveHyper/"

NAux = [3000]
alpha = 1
mode = "cluster"
clusters = 1
SimuIdentifier = "FullMonitor_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"
SimuIdentifierMean = "AvMembVolt_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

## Find and open the files
os.chdir(WorkingDirectory)
ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in np.sort(glob.glob(SimuIdentifier))[15:32]:
	J1 = float(re.search('J=(.+?)_', file1).group(1))
	sigma1 = float(re.search('sigma=(.+?)_', file1).group(1))	
	I01 = float(re.search('I0=(.+?)_', file1).group(1))

	Times = np.loadtxt(file1,skiprows=1,usecols=0,ndmin=1)
	NeuronVoltage = np.loadtxt(file1,skiprows=1,usecols=1,ndmin=1)
	NeuronCurrent = np.loadtxt(file1,skiprows=1,usecols=2,ndmin=1)

	NeuronVoltage -= np.nanmean(NeuronVoltage)
	fs = 1./(Times[1]-Times[0])*1000
	ft = np.fft.rfft(NeuronVoltage)
	freqs = np.fft.rfftfreq(len(NeuronVoltage), (Times[1]-Times[0])/1000) # Get frequency axis from the time axis
	mags = abs(ft) # We don't care about the phase information here

	Corr = signal.correlate(NeuronVoltage,NeuronVoltage)
	Corr /= np.max(Corr)
	Index = int(np.argmax(Corr))
	Corr = Corr[Index:]

	print(file1)
	plt.plot(Times,NeuronVoltage);plt.show();
	plt.plot(Corr);plt.show();quit()
	ResultsList = [J1, I01,  sigma1, mags]

	ListsOfLists.append(ResultsList)

DataMatrix = np.row_stack(ListsOfLists)
np.savetxt("/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes/DataMatrixAutoCorrGuillemHyper.txt",DataMatrix,fmt="%3.6f")

#J1, I01,  sigma1, SPC, R2, mean_Tnet = np.loadtxt("DataMatrixSPC.txt",skiprows=1,unpack=True)



