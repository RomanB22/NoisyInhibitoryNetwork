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

ChosenModel = 1 #Models = {0:'IzhiTypeI', 1:'IzhiTypeII', 2:'LIFModel', 3:'AdExModel', 4:'GuillemModel'}
ChosenSubmodelAdEx = 2 #SubModelsAdEx = {0:'TypeI', 1:'TypeIISaddleNode', 2:'TypeIIAndronovHopf'}
ThetaInput = False
FirstBisection = False # First run calculating Iext for N=800. This results serve as input for the second bisection search

VariableDrive=False
alphaAux = [1]# Scaling for the synaptic times
NAux = [3000]#[800,1000,1400,1800,2200,3000]
Nu0Aux = [30.]
ConnectionProbAux = [1]
InCond = ['cluster',1] #'random' or 'cluster' and second element is number of clusters
f = 8.*b2.Hz#8*b2.Hz # Frequency of the theta drive
IntegrationStep = 0.05 #in ms IntegrationStepGuillem = 0.01 and AdEx Type II Chattering (Andronov-Hopf Bifurcation)
transientTime = 1500*b2.ms#2000*b2.ms #Time until the network reaches the stationary state
simulationEnd = 1800*b2.ms#2300*b2.ms #Time of simulation end
PlotingTime = 300 #Time which is going to be plotted (in ms). Between 200 and 400 ms is okay. The interval goes between [simulationEnd-PlotingTime, simulationEnd]
SavingTime = 300*b2.ms #Time interval of the variables measured in the simulation to save. Be careful to use so much time here, because the size of all the generated .csv files could be huge. Maximum around 400 ms is okay. Interval is [simulationEnd-SavingTime, simulationEnd]
Neurons2Save = 1 #How many neurons to fully follow across the simulation time. Same as before, too much neurons creates huge .csv files
ChosenNeuronModelGuillem = 57 # MAximum membrane time constant, equal to 6.96 ms



#sigmaAux = [0.295065,0.334912,0.385459,0.400528,0.412271,0.472560,0.481536,0.563490,0.580941,0.649984, 0.665731,0.676694,0.804122,0.816118,0.899528,0.915436,1.092833,1.093559,1.131873,1.308772,1.315843, 1.555640,1.601277,1.828104,1.900892,2.050194,2.181402,2.429551,2.541174,2.804623,2.985417,3.066193, 3.292658,3.510757,3.693642,3.717038,3.820162,3.844092,3.923288,4.283818,5.037294,5.275650, 5.510934,5.785305,6.016972, 6.266655, 7.055545, 7.885936, 8.106726, 8.318635, 8.439412,8.810876]# IzhiType I 17Hz
#sigmaAux = [0.483636, 0.570125, 0.638976, 0.689093, 0.764012, 0.814777, 0.838117, 0.930955, 1.022948, 1.158297, 1.289844, 1.297850, 1.415989, 1.593778, 1.601389, 1.694786, 1.875540, 1.983402, 2.045293, 2.168014, 2.301241, 2.363165, 2.565969, 2.765484, 2.790878, 2.860648, 3.322227, 3.456152, 3.659212, 4.011020, 4.363336, 4.513535, 5.054940, 5.091141, 5.525228, 5.956527, 6.326659, 6.801862, 7.010915, 7.819147, 7.862861, 8.011350, 8.572892, 9.135147, 9.846626, 10.077171, 10.195416, 11.129005, 12.015768, 12.264251, 12.916860, 13.859298]# IzhiType I 30Hz
#sigmaAux = [0.271629, 0.291760, 0.296040, 0.335151, 0.470311, 0.563453, 0.568168, 0.799016, 0.919328, 1.122135, 1.398076, 1.400737, 1.760990, 2.114614, 2.171465, 2.529612, 2.989160, 3.047412, 3.562167, 3.584420, 3.811478, 4.099214, 4.271406, 4.882835, 5.010998, 5.130253, 5.761119, 5.958992, 6.014744, 6.511094, 7.010976, 7.010993, 7.355476, 7.784237, 8.344548, 8.451289, 8.531920, 8.918726, 9.371273, 9.526315, 9.785418, 9.801577, 10.642824, 11.099991, 11.126564, 11.519486, 11.665503, 12.033967, 12.409008, 12.779753, 13.318574, 13.628352]# Izhi Type II 17Hz
#sigmaAux = [0.338681, 0.400542, 0.454328, 0.556936, 0.563402, 0.585513, 0.804116, 0.804160, 0.949523, 1.308824, 1.555640, 1.579664, 1.878218, 2.394950, 2.548723, 2.855702, 3.216371, 3.572191, 4.040815, 4.226624, 4.589653, 5.183053, 5.260968, 5.803633, 5.956589, 6.672322, 7.197255, 7.356955, 7.581899, 8.464015, 8.760989, 8.843774, 9.581037, 10.132550, 11.022299, 11.048026, 11.439866, 11.920506, 12.847341, 13.341437, 13.456627, 14.420989, 15.118539, 15.813417, 15.976792, 16.094522, 16.192218, 16.754172, 18.057658, 18.117820, 18.216689, 18.777603]# Izhi Type II 30Hz
#sigmaAux = [0.070667, 0.085122, 0.086313, 0.099983, 0.114523, 0.136572, 0.138703, 0.161488, 0.192597, 0.266990, 0.277272, 0.278896, 0.331170, 0.339217, 0.394909, 0.460949, 0.480616, 0.546722, 0.562039, 0.768489, 0.803257, 0.945982, 1.068114, 1.087630, 1.104154, 1.126011, 1.192244, 1.240815, 1.364802, 1.400876, 1.622264, 1.862384, 1.881698, 1.941361, 2.145182, 2.318136, 2.344583, 2.358566, 2.644528, 2.811673, 3.178262, 3.225820, 3.307797, 3.674465, 3.797322, 3.970918, 4.199251, 4.332499, 4.445287, 4.591893, 4.810734, 5.022591] #AdEx Type II 17Hz
sigmaAux = [0.193379, 0.196648, 0.231697, 0.238631, 0.272892, 0.284768, 0.320308, 0.341051, 0.384073, 0.446735, 0.510994, 0.510995, 0.535392, 0.631035, 0.647837, 0.669053, 0.788820, 0.800822, 0.823624, 1.067356, 1.098269, 1.132062, 1.312023, 1.521582, 1.554428, 1.586002, 1.830308, 1.831646, 1.877379, 2.162752, 2.305446, 2.572358, 2.607653, 2.666924, 3.077807, 3.144101, 3.296498, 3.606819, 3.696050, 3.823642, 4.337391, 4.546868, 4.566507, 4.581596, 5.319661, 5.357002, 5.793489, 6.183745, 6.397163, 6.816094, 7.089126, 7.323615] #AdEx Type II 30Hz

JAux = [1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.2, 4.0, 5.0, 6.3, 7.9, 8.0, 10.0, 12.6, 13.0, 15.8, 18.0, 20.0, 23.0, 25.1, 28.0, 31.6, 33.0, 38.0, 39.8, 43.0, 48.0, 50.1, 53.0, 58.0, 63.0, 63.1, 68.0, 73.0, 78.0, 79.4, 83.0, 88.0, 93.0, 98.0, 100.0, 103.0, 108.0, 113.0, 118.0, 123.0, 128.0, 133.0, 138.0, 143.0, 148.0, 153.0]

ParametersListAux = [(sLoop, jLoop, aLoop, NLoop, CPLoop, Nu0Loop) for (sLoop,jLoop) in zip(sigmaAux,JAux) for aLoop in alphaAux for NLoop in NAux for CPLoop in ConnectionProbAux for Nu0Loop in Nu0Aux]

tau_m = 10.*b2.ms

# For Tigerfish use
maxArray = 1000 # Number of parallel array jobs sended to SLURM
RealizationNumber = int(sys.argv[2])

current_param = int(sys.argv[1])+RealizationNumber*maxArray

ParameterList = ParametersListAux[current_param]

WorkingDirectory = os.getcwd()

