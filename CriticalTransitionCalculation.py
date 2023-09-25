""" Critical Analysis around the transition point from the stable synchronous oscillation state to stable asynchronous state

In this program we read the Summary_ files for a determined set of initial conditions, alpha and connection probability, and create a 4D-array with all the values of [J, sigma, I0, N, SynchronyIndex]
In the asymptotic limit N -> infinite, the synchrony index Chi ~ Chi_inf + delta_x/sqrt(N) + O(1/N), with delta_x the first order finite-size correction.
Then, for each triplet [I0,J,sigma] I fit the function Chi(N) = Chi_inf + delta_x/sqrt(N), and obtain for each triplet [I0,J,sigma] the value of Chi_inf(J,sigma,I0)
For a fixed J and I0, near the transition Chi_inf behaves as A*(sigma_crit - sigma)^(1/2) for sigma < sigma_crit and 0 if sigma > sigma_crit
So, after all the process, I ended with the value of critical amount of noise for each J and I0 that destabilizes the asynchronous regime sigma_crit(J,I0)
Once the instability has been located, we compute the frequency of the population oscillations at the onset of the instability. To this end, we simulate the network for a noise level sigma ~ sigma_crit. A good estimate of the oscillation frequency is provided by looking at the autocorrelation of the population average of the membrane potentials of the neurons. If the bifurcation at sigma_crit is supercritical, the frequency estimates on the two sides of the transition are similar. If the bifurcation is subcritical, they may differ significantly. In that case, the frequency of the unstable mode has to be determined in the vicinity of the transition, but on the side where the asynchronous state is still stable

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
import os
import re

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

SavingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes" #os.getcwd
()
WorkingDirectory = os.path.dirname(os.path.normpath(SavingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = WorkingDirectory + CSVFolder
WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_IzhiTyIIVarDrive/"

alpha = 1
mode = "cluster"
clusters = 1
nu0 = 30
SimuIdentifier = "Summary_*J=1.6*-55.0.*"#"Summary_*_alpha=" + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"


sigmaAux = np.concatenate ( (np.around(np.linspace(0.011,10.011,51, endpoint=True),4) , np.around(np.linspace(10.011,15.011,11, endpoint=True),4)) )#np.around(np.linspace(0.011,10.011,51, endpoint=True),4)                                   
#JAux = [1.,20.,40.,60.,80.,100.,120.,140.,160.,180.,250.,300.]    
JAux =[1.6,79.4]# [1.6,79.4]
I0Aux= [0.,5.,10.,16.,20.,25.,30.,2.03125,3.59375,3.984375,4.2773,7.685547,15.29297]  

## Find and open the files
os.chdir(WorkingDirectory)
ListsOfLists = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))
#print(NumberOfSimulations);quit()

for file1 in glob.glob(SimuIdentifier):
    with open(file1,'r') as f:
        FullList = []
        for key,group in it.groupby(f,lambda line: line.startswith('Simulation Summary') or \
        line.startswith(' alpha') or line.startswith(' Connection') or\
        line.startswith(' Avg.') or line.startswith(' InitialConditions') or line.startswith(' Binder Cumulant')):
            if not key:
                List = []
                for line in group:
                    Number = re.findall(r'\d+.\d+', line) # find number of digits through regular expression
                    if np.size(Number)==0: Number = [0]
                    List.append(float(Number[0]))
                FullList.append(List)
        ResultsList = list(it.chain.from_iterable(FullList))# ResultsList = [J, sigma, N, Chi]
        ListsOfLists.append(ResultsList)

"""
ListsOfLists has the four values of [I0,J, sigma, N, Chi] for each one of the simulations
So ListsOfLists[0] shows the results for the first simulation
   ListsOfLists[indexSimulation][0] = I0 for simulation = IndexSimulation
   ListsOfLists[indexSimulation][1] = J for simulation = IndexSimulation
   ListsOfLists[indexSimulation][2] = sigma for simulation = IndexSimulation
   ListsOfLists[indexSimulation][3] = N for simulation = IndexSimulation
   ListsOfLists[indexSimulation][4] = Chi for simulation = IndexSimulation
"""
DataMatrix = np.row_stack(ListsOfLists)

if ConductanceBased == True:
    DataMatrix[:,1] = DataMatrix[:,1]/1e3

print(DataMatrix);quit()
plotFlag = False # Whether to plot or not the synchrony index versus population size
plotSigmaCritFlag = False # Whether to plot or not the Chi_inf versus sigma

def AsymptoticSynchrony(N,Chi_inf,Delta_x):
    return Chi_inf + Delta_x/np.sqrt(N)
    
def CriticalBehavior(x,A,sigma_crit): # B*(sigma_crit-x)**1.5
    return np.where(x < sigma_crit, A*(sigma_crit-x)**0.5, 0)
    
def objective(params, x, y):
    return np.sum(np.abs(CriticalBehavior(x, *params) - y))

NumberOfI0s = np.shape(I0Aux)[0]
NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]

#print(I0Aux,JAux)
#HeatMapChi_inf = -1.*np.ones((NumberOfSigmas,NumberOfJs))
#indexJ = 0


for I0 in I0Aux:
    CriticalLine = []
    if I0<0:
        IndexI0 = np.argwhere(abs(DataMatrix[:,0]+I0)<1e-2)
    else:
        IndexI0 = np.argwhere(abs(DataMatrix[:,0]-I0)<1e-2)
    for J in JAux:
        indexSigma = 0
        Chi_inf_sigma = [] # For each J, (sigma, Chi_inf, Error in Chi_inf) [I0,J, sigma, N, Chi]
        for sigma in sigmaAux:
            IndexJ = np.argwhere(abs(DataMatrix[:,1]-J)<1e-2)
            IndexSigma = np.argwhere(abs(DataMatrix[:,2]-sigma)<1e-4)
            Indexes = np.intersect1d(IndexI0,np.intersect1d(IndexJ,IndexSigma))
            #print(Indexes,np.intersect1d(IndexJ,IndexSigma),IndexJ,IndexSigma,IndexI0)

            if Indexes.size is not 0:
                NVector = DataMatrix[Indexes,3]
                Chi_N = DataMatrix[Indexes,4]
                if(sum(np.shape(Chi_N))>2): ## This part is only for testing. When all the simulations where done this condition is always true
                    popt, pcov = curve_fit(AsymptoticSynchrony, NVector, Chi_N, bounds = ([0.,-np.inf],[1,np.inf]))
                    Chi_inf_sigma.append([sigma,popt[0],pcov[0,0]])
                    if plotFlag:
                        plt.figure()
                        plt.plot(NVector,Chi_N,'b*',label='simulation')                
                        AuxVector = np.linspace(800,3000,400)
                        plt.plot(AuxVector, AsymptoticSynchrony(AuxVector, *popt), 'r-', label='fit: $\chi_\infty$=%5.3f, $\delta_x$=%5.3f' % tuple(popt))
                        plt.legend()
                        plt.title("J = %d, $\sigma$ = %.2f" % (J,sigma))
                        plt.xlabel('Network size')
                        plt.ylabel('$\chi(N)$')
                        NameFig = SavingDirectory + "/ChiN/" + "SynchrIndx_J=%d_sigma=%.2f" % (J,sigma)
                        plt.savefig(NameFig+".png")
                        plt.savefig(NameFig+".svg")
                        plt.show()
                        plt.close()
                    DataChi_inf = np.row_stack(Chi_inf_sigma)
                    #HeatMapChi_inf[indexSigma,indexJ] = popt[0]
                    #print(popt[0])
            #indexSigma += 1
        #indexJ += 1
        
        if ChosenModel==0 or ChosenModel==2 or ChosenModel==3 or ChosenModel==5 or ChosenModel==6:
            guess = [1.1,2.1]
            VectorX = DataChi_inf[:30,0]
            VectorY = DataChi_inf[:30,1]
            res = minimize(objective, guess, args=(VectorX,VectorY), method='Nelder-Mead')
            if UniformReversal==False and abs(Esyn/b2.mV - -75.)<1e-3:
                if J > 70 and J<180:
                    VectorX = DataChi_inf[20:,0]
                    VectorY = DataChi_inf[20:,1]
                    guess = [1.1,5.1]
                elif  J >= 180 and J<250:
                    VectorX = DataChi_inf[25:,0]
                    VectorY = DataChi_inf[25:,1]
                    guess = [1.1,6.1]
                elif  J >= 250:
                    VectorX = DataChi_inf[25:,0]
                    VectorY = DataChi_inf[25:,1]
                    guess = [1.1,10.1]                
                res = minimize(objective, guess, args=(VectorX,VectorY), method='Nelder-Mead')    

        else:
            #index = np.argmin(DataChi_inf[1:,1])
            #indexBegin = max(0,index-15)
            #indexEnd = min(index+10,np.shape(DataChi_inf)[0])
            #print(J,index,indexBegin,indexEnd)
            guess = [1.1,3.1]
            res = minimize(objective, guess, args=(DataChi_inf[:,0],DataChi_inf[:,1]), method='Nelder-Mead')
            if I0<10 and J<40: res = minimize(objective, [0.01,4], args=(DataChi_inf[:,0],DataChi_inf[:,1]), method='Nelder-Mead')    
            if I0<10 and J==40: res = minimize(objective, [0.01,4], args=(DataChi_inf[:,0],DataChi_inf[:,1]), method='Nelder-Mead')    
            if J==1: res = minimize(objective, [0.01,0.5], args=(DataChi_inf[1:15,0],DataChi_inf[1:15,1]), method='Nelder-Mead')
            if J==20 and abs(I0-75.)<1e-1: res = minimize(objective, [0.7,3.8], args=(DataChi_inf[:,0],DataChi_inf[:,1]), method='Nelder-Mead')

        if plotSigmaCritFlag:
            plt.figure()
            plt.plot(DataChi_inf[1:,0],DataChi_inf[1:,1],'bx',label='data')
            x = np.linspace(DataChi_inf[0,0],max(DataChi_inf[0:,0]),400)
            plt.plot(x,CriticalBehavior(x,*(res.x)),'r-',label='fit: $\sigma_c$ = %1.3f' % res.x[1])
            plt.title("J = %d" % J)
            plt.xlabel('$\sigma$ (mV)')
            plt.ylabel('$\chi_\infty$')
            #plt.xscale('log')
            #plt.yscale('log')
            #plt.xlim([0.01,10.011])
                    #plt.ylim([0,1])
            plt.legend()
            NameFig = SavingDirectory + "/ChiInf/" + "{Folder}ChiInf_J={J}_I0={I0}".format(Folder=CSVFolder[1:-1],J=J,I0=I0)
            plt.savefig(NameFig+".png")
            plt.savefig(NameFig+".svg")
            plt.show()
            plt.close()

        CriticalLine.append([I0,res.x[1],J])
    MatrixCritLine = np.row_stack(CriticalLine)
    np.savetxt(SavingDirectory + "/16CriticalTransitionLine"+CSVFolder[1:-1]+"_I0"+str(I0)+"Hz.csv",MatrixCritLine,delimiter = " ", fmt="%3.1f %2.6f %3.1f")



