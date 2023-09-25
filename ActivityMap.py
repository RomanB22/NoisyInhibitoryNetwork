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

sigmaAux = np.concatenate ( (np.around(np.linspace(0.011,10.011,51, endpoint=True),4) , np.around(np.linspace(10.011,15.011,11, endpoint=True),4)) )#np.around(np.linspace(0.011,10.011,51, endpoint=True),4)                                   JAux = [1.,20.,40.,60.,80.,100.,120.,140.,160.,180.,250.,300.]                                                          I0Aux= [30.,40.,45.,50.,55.,60.,65.,70.,75.,80.]                                                                                                                                                                                                
JAux =[79.4]# [1.6,79.4]
I0Aux= [0.,5.,10.,16.,20.,25.,30.,2.03125,3.59375,3.984375,4.2773,7.685547,15.29297]   

SavingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes" #os.getcwd()
WorkingDirectory = os.path.dirname(os.path.normpath(SavingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = WorkingDirectory + CSVFolder
WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_IzhiTyIIVarDrive/"

print(CSVFolder,SavingDirectory)

NAux = [3000]
alpha = 1
mode = "cluster"
clusters = 1
nu0 = 17
#SimuIdentifier = "Raster_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"
#SimuIdentifier = "Raster_*J=1.6*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"
SimuIdentifier = "Raster_*J=79.4*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"

def elements(array):
	return array.ndim and array.size

## Find and open the files
os.chdir(WorkingDirectory) #Change working directory

ListsOfLists = []
NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
    print(file1)
    J1 = re.search('J=(.+?)_', file1).group(1)
    sigma1 = re.search('sigma=(.+?)_', file1).group(1)
    I1 = re.search('I0=(.+?)_', file1).group(1)
    Raster =  np.loadtxt(file1,skiprows=1,usecols=1)
    print(np.loadtxt(file1) )
    ID = np.unique(Raster)
    Results = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float),np.array([I1]).astype(np.float),np.shape(ID)[0]]
    #if np.shape(ID)[0]==NAux[0]:
    #	Results = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float),1]
    #else:
    #	Results = [np.array([J1]).astype(np.float), np.array([sigma1]).astype(np.float),0]
    ListsOfLists.append(Results)

DataMatrix = np.row_stack(ListsOfLists)

NumberOfJs = np.shape(JAux)[0]
NumberOfSigmas = np.shape(sigmaAux)[0]
HeatMapActivity = -1.*np.ones((NumberOfSigmas,NumberOfJs))

indexJ = 0
for J in JAux:
    indexSigma = 0
    for sigma in sigmaAux:
        IndexJ = np.argwhere(abs(DataMatrix[:,0]-J)<0.1)
        IndexSigma = np.argwhere(abs(DataMatrix[:,1]-sigma)<0.0001)
        Indexes = np.intersect1d(IndexJ,IndexSigma)
        if Indexes.size is not 0:
            HeatMapActivity[indexSigma,indexJ] = DataMatrix[Indexes,3][0]
        indexSigma += 1
    indexJ += 1
dataframeActivity = pandas.DataFrame(HeatMapActivity,index=sigmaAux,columns=JAux)
print(dataframeActivity)
#dataframeActivity.to_csv(SavingDirectory +"/ActivityMask_16_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1])) 
dataframeActivity.to_csv(SavingDirectory +"/ActivityMask_794{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1])) 
#dataframeActivity.to_csv(SavingDirectory +"/ActivityMask_{nu}_{Folder}.csv".format(nu=nu0,Folder=CSVFolder[1:-1]))
