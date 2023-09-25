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

SavingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/ScriptsAnalysis/CriticalTransition/InputPlanes" #os.getcwd()
WorkingDirectory = os.path.dirname(os.path.normpath(SavingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = os.path.dirname(os.path.normpath(WorkingDirectory))
WorkingDirectory = WorkingDirectory + CSVFolder
WorkingDirectory = "/mnt/beegfs/home/rbarav/FullNetworkSimulations/NewPackageVersion/CSV_IzhiTyIIVarDrive/"

sigmaAux = np.concatenate ( (np.around(np.linspace(0.011,10.011,51, endpoint=True),4) , np.around(np.linspace(10.011,15.011,11, endpoint=True),4)) )#np.around(np.linspace(0.011,10.011,51, endpoint=True),4)
JAux = [1.,20.,40.,60.,80.,100.,120.,140.,160.,180.,250.,300.]
I0Aux= [30.,40.,45.,50.,55.,60.,65.,70.,75.,80.]
JAux =[1.6]# [1.6,79.4]
I0Aux= [0.,5.,10.,16.,20.,25.,30.,2.03125,3.59375,3.984375,4.2773,7.685547,15.29297]     

NAux = [3000]
alpha = 1
mode = "cluster"
clusters = 1
nu0 = 17
SimuIdentifier = "Summary_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

#SimuIdentifier = "Summary_*J=1.6*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"
SimuIdentifier = "Summary_*J=1.6*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"


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

##############
SimuIdentifier = "FullMonitor_*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"
#SimuIdentifier = "FullMonitor_*J=1.6*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"
SimuIdentifier = "FullMonitor_*J=1.6*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"



ListsOfLists = []

PermanenceTimeInEachPhaseSpace = []

NumberOfSimulations = len(glob.glob(SimuIdentifier))

for file1 in glob.glob(SimuIdentifier):
    J1 = re.search('J=(.+?)_', file1).group(1)
    sigma1 = re.search('sigma=(.+?)_', file1).group(1)
    Full = np.loadtxt(file1,skiprows=0)
    Text = "Summary_*J=79.4*N={N}_alpha=".format(N=NAux[0]) + str(alpha) + '*_InitCond=' + str(mode) + str(clusters) + "*-55.0.csv"#"Summary_J={J}*_sigma={sigma}_*N={N}_alpha=".format(J=J1,sigma=sigma1,N=NAux[0]) + str(alpha) + '*_nu0=' + str(nu0) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"
    I0Map = glob.glob(Text)[0]
    with open(I0Map,'r') as f:
        FullList = []
        for key,group in it.groupby(f,lambda line: line.startswith(' I0 =')):
            if key:
                for line in group:
                    Number = np.array(re.findall(r'\d+.\d+', line)).astype(np.float) # find number of digits through regular expression
                    if np.size(Number) == 0:
                        Number = [0.]
    print(file1)
    if np.shape(Full)[1]==16:    
        Irec = np.loadtxt(file1,skiprows=0,usecols=11)# Start at index 11 to 15 (for N different of 800). For N=800 starts at index 6 to 10
    if np.shape(Full)[1]==11:
        Irec = np.loadtxt(file1,skiprows=0,usecols=6)
    if np.shape(Full)[1]==4:
        Irec = np.loadtxt(file1,skiprows=0,usecols=3)
    ITotDiscretized = Irec+Number[0]

    #ITotDiscretized[ITotDiscretized<1.795]=1
    #ITotDiscretized[np.logical_and(ITotDiscretized>=1.795, ITotDiscretized<2.625)]=2
    #ITotDiscretized[np.logical_and(ITotDiscretized>=2.625,  ITotDiscretized<4.225)]=3
    #ITotDiscretized[ITotDiscretized>4.225]=4

    #ITotDiscretized[ITotDiscretized<0.129321*10/0.09]=1
    #ITotDiscretized[ITotDiscretized>=0.129321*10/0.09]=2

    #ITotDiscretized[ITotDiscretized<1.795]=1
    #ITotDiscretized[np.logical_and(ITotDiscretized>=1.795, ITotDiscretized<2.625)]=2
    #ITotDiscretized[np.logical_and(ITotDiscretized>=2.625,  ITotDiscretized<4.225)]=3
    #ITotDiscretized[ITotDiscretized>4.225]=4

    hist,bins = np.histogram(ITotDiscretized, bins=2,range=(0.9,2.1))
    IrecMean = np.mean(Irec)
    IrecMax = max(Irec)
    IrecMin = min(Irec)
    ResultsList = [np.array(J1).astype(np.float), np.array(sigma1).astype(np.float), IrecMean,IrecMax,IrecMin]
    ResultPermanence = [np.array(J1).astype(np.float), np.array(sigma1).astype(np.float), hist[0],hist[1]]
    ListsOfLists.append(ResultsList)
    PermanenceTimeInEachPhaseSpace.append(ResultPermanence)
DataMatrix2 = np.row_stack(ListsOfLists)
DataMatrixPermanence = np.row_stack(PermanenceTimeInEachPhaseSpace)

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

        IndexJ = np.argwhere(abs(DataMatrix2[:,0]-J)<0.1)
        IndexSigma = np.argwhere(abs(DataMatrix2[:,1]-sigma)<0.0001)
        Indexes = np.intersect1d(IndexJ,IndexSigma)

        if Indexes.size is not 0:
            HeatMapPermanenceSubSNP[indexSigma,indexJ] = DataMatrixPermanence[Indexes,2][0]
            HeatMapPermanenceBist[indexSigma,indexJ] = DataMatrixPermanence[Indexes,3][0]
            #HeatMapPermanenceSupraHopf[indexSigma,indexJ] = DataMatrixPermanence[Indexes,4][0]
            #HeatMapPermanenceSupraSN[indexSigma,indexJ] = DataMatrixPermanence[Indexes,5][0]
            Irec = DataMatrix2[Indexes,2][0]
            HeatMapIrecMean[indexSigma,indexJ] = Irec
            IrecMax = DataMatrix2[Indexes,3][0]
            HeatMapIrecMax[indexSigma,indexJ] = IrecMax
            IrecMin = DataMatrix2[Indexes,4][0]
            HeatMapIrecMin[indexSigma,indexJ] = IrecMin
        #print(J,sigma,HeatMapI0[indexSigma,indexJ]+HeatMapIrec[indexSigma,indexJ])
        indexSigma += 1
    indexJ += 1

dataframe = pandas.DataFrame(HeatMapI0+HeatMapIrecMean,index=sigmaAux,columns=JAux)
mask = dataframe == -300.
dataframe[mask] = np.nan
dataframe.index = dataframe.index.to_series().apply(lambda x: np.round(x,4))

dataframe2 = pandas.DataFrame(HeatMapI0+HeatMapIrecMax,index=sigmaAux,columns=JAux)
mask = dataframe2 == -300.
dataframe2[mask] = np.nan
dataframe2.index = dataframe2.index.to_series().apply(lambda x: np.round(x,4))


dataframeSubSNP = pandas.DataFrame(HeatMapPermanenceSubSNP,index=sigmaAux,columns=JAux)
dataframeBist = pandas.DataFrame(HeatMapPermanenceBist,index=sigmaAux,columns=JAux)
#dataframeSupraHopf = pandas.DataFrame(HeatMapPermanenceSupraHopf,index=sigmaAux,columns=JAux)
#dataframeSupraSN = pandas.DataFrame(HeatMapPermanenceSupraSN,index=sigmaAux,columns=JAux)

#dataframeSubSNP.to_csv(SavingDirectory+"/MapSubSNCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
#dataframeBist.to_csv(SavingDirectory+"/MapSupraSNCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))

#dataframeSubSNP.to_csv(SavingDirectory+"/MapSubSNPCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
#dataframeBist.to_csv(SavingDirectory+"/MapBistCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
#dataframeSupraHopf.to_csv(SavingDirectory+"/MapSupraHopfCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
#dataframeSupraSN.to_csv(SavingDirectory+"/MapSupraSNCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))

#print("Early end")
#quit()

if ChosenModel == 0:
    dataframe.to_csv(SavingDirectory+"/MapMeanInputCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
    df1 = dataframe >= 0.129321*10/0.09
    df2 = dataframe < 0.129321*10/0.09
    dataframe[df1] = 2
    dataframe[df2] = 1

    plt.figure()
    ax = sns.heatmap(dataframe.T.iloc[::-1],cmap="YlGnBu",square=False, vmin = 0.9, vmax = 2.1)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([1, 2])
    cbar.set_ticklabels(['Below SN', 'Above SN'])
    plt.clf()
        
    dataframe2.to_csv(SavingDirectory+"/MapMaxInputCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
    df1 = dataframe2 >= 0.129321*10/0.09
    df2 = dataframe2 < 0.129321*10/0.09
    dataframe2[df1] = 2
    dataframe2[df2] = 1

elif ChosenModel == 1:
    dataframe.to_csv(SavingDirectory+"/MapMeanInputCurrent16_{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
    df1 = dataframe > 0.4225*10.
    df2 = dataframe < 0.4225*10.
    df3 = dataframe > 0.2625*10.
    df4 = dataframe < 0.2625*10.
    df5 = dataframe > 0.1795*10.
    df6 = dataframe < 0.1795*10.
    dataframe2.to_csv(SavingDirectory+"/MapMaxInputCurrent16_{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))


    dataframe[df1] = 4
    dataframe[df2.mul(df3)] = 3
    dataframe[df4.mul(df5)] = 2
    dataframe[df6] = 1
    plt.figure()
    ax = sns.heatmap(dataframe.T.iloc[::-1],cmap="YlGnBu",square=False, vmin = 0.9, vmax = 4.1)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([1, 2, 3, 4])
    cbar.set_ticklabels(['Below SNP', 'Between SNP and Hopf', 'Between Hopf and SN', 'Above SN'])

elif ChosenModel==3:
    factor=10./np.array(CAdEx)
    if ChosenSubmodelAdEx==2: # AdEx AndronovHopf
        print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        dataframe.to_csv(SavingDirectory+"/MapMeanInputCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
        df1 = dataframe > factor*185.2901
        df2 = dataframe < factor*185.2901
        df3 = dataframe > factor*184.8613
        df4 = dataframe < factor*184.8613
        df5 = dataframe > factor*165
        df6 = dataframe < factor*165
        dataframe2.to_csv(SavingDirectory+"/MapMaxInputCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))


        dataframe[df1] = 4
        dataframe[df2.mul(df3)] = 3
        dataframe[df4.mul(df5)] = 2
        dataframe[df6] = 1

        plt.figure()
        ax = sns.heatmap(dataframe.T.iloc[::-1],cmap="YlGnBu",square=False, vmin = 0.9, vmax = 4.1)
        cbar = ax.collections[0].colorbar
        cbar.set_ticks([1, 2, 3, 4])
        cbar.set_ticklabels(['Below SNP', 'Between SNP and Hopf', 'Between Hopf and SN', 'Above SN'])

elif ChosenModel == 2:
        dataframe.to_csv(SavingDirectory+"/MapMeanInputCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
        df1 = dataframe >= 20
        df2 = dataframe < 20
        dataframe[df1] = 2
        dataframe[df2] = 1

        plt.figure()
        ax = sns.heatmap(dataframe.T.iloc[::-1],cmap="YlGnBu",square=False, vmin = 0.9, vmax = 2.1)
        cbar = ax.collections[0].colorbar
        cbar.set_ticks([1, 2])
        cbar.set_ticklabels(['Below SN', 'Above SN'])
        plt.clf()

        dataframe2.to_csv(SavingDirectory+"/MapMaxInputCurrent{N}_{nu}_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))
        df1 = dataframe2 >= 20
        df2 = dataframe2 < 20
        dataframe2[df1] = 2
        dataframe2[df2] = 1


plt.xlabel('$\sigma$ (mV)')
plt.ylabel('J (mV)')
plt.legend()
plt.tight_layout()
plt.title('$mean(I_0 + I_{syn})$')

#plt.savefig(SavingDirectory+"/HeatMaps/HeatMapMeanInputCurrent{N}_{nu}_{Folder}.svg".format(N=NAux[0],nu=nu0,Folder=CSVFolder[1:-1]))

