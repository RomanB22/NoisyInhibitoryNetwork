"""File for defining the simulation parameters


Parameters
----------

ParameterList : List of lists of parameters
   
	Each list has 6 elements: 
		ParameterList[0]: Noise strength (sigma), in mV. Watch out for the different normalization to each neuron model (LIF or Izhikevich) in order to make a fair comparison between models
		ParameterList[1]: Coupling strength (J), in mV
		ParameterList[2]: Synaptic times factor (alpha)
		ParameterList[3]: Number of neurons in the network (N)
		ParameterList[4]: Connection probability (epsilon). 1 means fully connected -without autapses-; 0 means unconnected. The connectivity pattern is random (Erdos-Renyi)
		ParameterList[5]: Desired average frequency (Nu_0). It is used in the I0 bisection search, and to keep a track in the simulation Summary file

maxArray : Number of parallel array jobs sended to SLURM. Should coincide with the array size in the .slurm file

transientTime : Time until the network reaches the stationary state

simulationEnd : Time of simulation end

PlotingTime : Time which is going to be plotted (in ms). Between 200 and 400 ms is okay. The interval goes between [simulationEnd-PlotingTime, simulationEnd]

SavingTime : Time interval of the variables measured in the simulation to save. Be careful to use so much time here, because the size of all the generated .csv files could be huge. Maximum around 400 ms is okay. Interval is [simulationEnd-SavingTime, simulationEnd]

Neurons2Save : How many neurons to fully follow across the simulation time. Same as before, too much neurons creates huge .csv files

FiguresFolder : Folder where save plots. 

CSVFolder : Folder where save the .csv files. 

I0Folder : Folder where save and load the .csv files with the results for the DC input 

f : Frequency of the theta drive

"""

b2.seed(7896) # Set random seed to have reproducible results!

ChosenModel = 4 #Models = {0:'IzhiTypeI', 1:'IzhiTypeII', 2:'LIFModel', 3:'AdExModel', 4:'GuillemModel'}
ChosenSubmodelAdEx = 2 #SubModelsAdEx = {0:'TypeI', 1:'TypeIISaddleNode', 2:'TypeIIAndronovHopf'}
ThetaInput = True
FirstBisection = False # First run calculating Iext for N=800. This results serve as input for the second bisection search


alphaAux = [1]# Scaling for the synaptic times
NAux = [800]#[800,1000,1400,1800,2200,3000]
Nu0Aux = [17.]
ConnectionProbAux = [1]
InCond = ['cluster',1] #'random' or 'cluster' and second element is number of clusters
f = 8.*b2.Hz#8*b2.Hz # Frequency of the theta drive
IntegrationStep = 0.05 #in ms IntegrationStepGuillem = 0.01 and AdEx Type II Chattering (Andronov-Hopf Bifurcation)
transientTime = 1200*b2.ms#1500*b2.ms #900#Time until the network reaches the stationary state
simulationEnd = 3100*b2.ms#1800*b2.ms #1200#Time of simulation end
PlotingTime = 1900 #Time which is going to be plotted (in ms). Between 200 and 400 ms is okay. The interval goes between [simulationEnd-PlotingTime, simulationEnd]
SavingTime = 1900*b2.ms #Time interval of the variables measured in the simulation to save. Be careful to use so much time here, because the size of all the generated .csv files could be huge. Maximum around 400 ms is okay. Interval is [simulationEnd-SavingTime, simulationEnd]
Neurons2Save = 5 #How many neurons to fully follow across the simulation time. Same as before, too much neurons creates huge .csv files
ChosenNeuronModelGuillem = 57 # MAximum membrane time constant, equal to 6.96 ms

sigmaAux = [0.1334,0.1884,0.5110,1.2589,3.7610,5.7610,9.7610]#np.sort(np.concatenate((np.around(np.logspace(-2,1,41,endpoint=True),4),np.linspace(0.011,10.011,41, endpoint=True))))

JAux = [10,20,100,153]#np.sort(np.concatenate((np.around(np.logspace(0,2,21,endpoint=True),1),np.linspace(3,153,31, endpoint=True))))

ParametersListAux = [(sLoop, jLoop, aLoop, NLoop, CPLoop, Nu0Loop) for sLoop in sigmaAux for jLoop in JAux for aLoop in alphaAux for NLoop in NAux for CPLoop in ConnectionProbAux for Nu0Loop in Nu0Aux]

tau_m = 10.*b2.ms

# For Tigerfish use
maxArray = 1000 # Number of parallel array jobs sended to SLURM
RealizationNumber = int(sys.argv[2])

current_param = int(sys.argv[1])+RealizationNumber*maxArray

ParameterList = ParametersListAux[current_param]

WorkingDirectory = os.getcwd()

