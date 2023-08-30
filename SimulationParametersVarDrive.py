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
import brian2 as b2
import sys, os
import numpy as np


b2.seed(7896) # Set random seed to have reproducible results!

ChosenModel = 1 #Models = {0:'IzhiTypeI', 1:'IzhiTypeII', 2:'LIFModel', 3:'AdExModel', 4:'GuillemModel', 5:'ModifiedIzhiTyII', 6:'WangBuzsaki'}
ChosenSubmodelAdEx = 2 #SubModelsAdEx = {0:'TypeI', 1:'TypeIISaddleNode', 2:'TypeIIAndronovHopf'}
ThetaInput = False
FirstBisection = False # First run calculating Iext for N=800. This results serve as input for the second bisection search
VariableDrive = True
ConductanceBased = False; 
UniformReversal = False;
Esyn = -55.*b2.mV#-75.*b2.mV; #shunting -55 mV, hyperpolarizing = -75 mV

alphaAux = [1]# Scaling for the synaptic times
NAux = [3000]#[800,1400,2200,3000]
Nu0Aux = [30.]
ConnectionProbAux = [1]
InCond = ['cluster',1]#'random' or 'cluster' and second element is number of clusters
f = 8*b2.Hz # Frequency of the theta drive
IntegrationStep = 0.05 #in ms IntegrationStepGuillem = 0.01 and AdEx Type II Chattering (Andronov-Hopf Bifurcation)
transientTime = 1100*b2.ms#1500*b2.ms#1600*b2.ms #Time until the network reaches the stationary state
simulationEnd = 1500*b2.ms#1900*b2.ms#2200*b2.ms #Time of simulation end
PlotingTime = 200 #Time which is going to be plotted (in ms). Between 200 and 400 ms is okay. The interval goes between [simulationEnd-PlotingTime, simulationEnd]
SavingTime = 400*b2.ms #Time interval of the variables measured in the simulation to save. Be careful to use so much time here, because the size of all the generated .csv files could be huge. Maximum around 400 ms is okay. Interval is [simulationEnd-SavingTime, simulationEnd]
Neurons2Save = 1 #How many neurons to fully follow across the simulation time. Same as before, too much neurons creates huge .csv files
ChosenNeuronModelGuillem = 57 # 43 is Ananth using

tau_m = 10.*b2.ms

#sigmaAux = np.around(np.linspace(0.011,10.011,51, endpoint=True),4)
#JAux = np.around(np.logspace(0,2,21,endpoint=True),1)#Plane I0xJ
#JAux = [1.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.] #[1.,20.,40.,60.,80.,100.,120.,140.]#np.around(np.logspace(0,2,11,endpoint=True),1)#Plane I0xsigma
#I0Aux = [0.,5.,10.,15.,20.,25.,30.,35.,40.]

#Conductance based
sigmaAux = np.around(np.linspace(10.011,15.011,11,endpoint=True))
#np.sort( np.concatenate( (np.around(np.linspace(10.011,15.011,11,endpoint=True),4),np.around(np.linspace(0.011,10.011,51, endpoint=True),4)) ) )
JAux = [100.]#[1.,20.,40.,60.,80.,100.,120.,140.,160.,180.,250.,300.]
I0Aux= [1.,4.,7.,10.,13.,16.,19.,22.,25.,28.]#[30.,40.,45.,50.,55.,60.,65.,70.,75.,80.]

ParametersListAux = [(sLoop, jLoop, aLoop, NLoop, CPLoop, I0Loop) for sLoop in sigmaAux for jLoop in JAux for aLoop in alphaAux for NLoop in NAux for CPLoop in ConnectionProbAux for I0Loop in I0Aux]


# For Tigerfish use
maxArray = 1000 # Number of parallel array jobs sended to SLURM
RealizationNumber = int(sys.argv[2])

current_param = int(sys.argv[1])+RealizationNumber*maxArray

ParameterList = ParametersListAux[current_param]

WorkingDirectory = os.getcwd()

