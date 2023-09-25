"""Main file for the simulation of the Izhikevich type 2 resonator network with fully connected geometry, according to the model of Brunel & Hansel 2006

Simulation of the Izhikevich network with fully connected geometry, according to the model of Brunel & Hansel 2006. The main goal of the simulation is to understand the influence of the connectivity strength and the noise in the stability of the asynchronous state, and the transitions between the different synchronous states and the asynchronous state.
There are two main instabilities of the asynchoronous state: 'clustering' instabilities (predominant in the coupled oscillator regions -low noise and low coupling strength-) and the 'oscillatory rate instability'/'stochastic population oscillator' (predominant in the high noise and high coupling region). 

Author 
------
Roman Baravalle <romanbaravalle@gmail.com>

Version
-------
0.2

Date
----
10-01-2021

"""
execfile("ImportingPackages.py")
execfile("SimulationParametersMain.py")
execfile("AuxiliarFunctions.py")
execfile("NeuronModels.py")
execfile("SetupModels.py")

## Initial time of the program
print("Initializing program...")
if ChosenModel==3:
	print("Model: " + Models[ChosenModel]+" "+SubModelsAdEx[ChosenSubmodelAdEx])
else:
	print("Model: " + Models[ChosenModel])
t1 = time.time() # Storing the begining time of the simulation

"""
Creating the network
Here we create the network for a set of parameters. We want to adjust the external DC input I0 in a way that the average firing rate of the network is around a desired fixed value

"""

b2.start_scope() # This is here for running separate simulations in the same notebook
b2.defaultclock.dt = IntegrationStep * b2.ms # This is the timestep Brunel uses in his paper

"""Assign values to each variable of the network"""
Nu_0 = ParameterList[5]*b2.Hz # desired average population frequency
DeltaNu_0 = 0.2*b2.Hz # expected error in the average mean firing rate (only useful for the SolveI0 routine)
J = ParameterList[1]*b2.mV
sigma = ParameterList[0]*b2.mV
N=ParameterList[3]

"""Create NeuronGroup"""
print("Creating NeuronGroup...")

if ChosenModel==4:
	Network = b2.NeuronGroup(N, GuillemModel, threshold=threshold_cond, refractory=.7*b2.ms, method="euler") # be careful with spike counting: should consider one spike / spike and not many when above v_thr. thats why we add refractoriness
else:
	Network = b2.NeuronGroup(N, Model, threshold=threshold_cond, reset=reset_cond, method="euler")

"""Set initial conditions"""
mode = InCond[0]
clusters = InCond[1]
if ChosenModel==4:
	Network.EL, Network.gL, Network.tauv, Network.Cap = ELs, gLs, tauvs, Caps;
	Network.gNa, Network.gKv3, Network.thm1, Network.thh2, Network.thn1 = gNas, gKv3s, thm1s, thh2s, thn1s;
	Network.gKv1, Network.tha1 = gKv1s, tha1s;
	"""
	# Set initial conditions for v and gating variables
	Network.v, Network.m, Network.h, Network.n = ELs+3.*np.random.randn()*b2.mV, .01+.1*np.random.rand(), .99-.1*np.random.rand(), .01+.1*np.random.rand(); # Random initial values for vs and gating variables
	Network.a = .01+.1*np.random.rand();
	"""
	Network.v, Network.m, Network.h, Network.n = ELs, .01, .99, .01; # Constant initial values for vs and gating variables
	Network.a = .01;

else:
	v0, u0 = SetInitialVoltages(N, mode, clusters, v_th=v_th, v_r=v_r, v0=v0, u0=u0)
	Network.v = v0*b2.mV

	if ChosenModel==0 or ChosenModel==1 or ChosenModel==5:
		Network.u = u0*uUnits
	elif ChosenModel==3:
		Network.w = u0*wUnits

"""
Set synapses

Parameters
----------
tau_l : latency time of the postsynaptic current, in ms
alpha : relation between synaptic times

"""
alpha = ParameterList[2] 
tau_d = alpha*6*b2.ms
tau_l = alpha*b2.ms
tau_r = alpha*b2.ms
S = b2.Synapses(Network,Network,model=synaps_eqs,on_pre=synaps_action,method="euler")

"""All-to-all connectivity"""
epsilon = ParameterList[4] # p=1 means 100% probability of connection
S.connect(p=epsilon,condition="i!=j")

S.delay = tau_l # Set delay
S.DeltaX = tau_m/tau_r # Set postsynaptic increment

t2 = time.time() # Calculating initialization time
print("\tInitialization took %.3f sec" % (t2-t1))

"""Store initial state for restoring during bisection search and full simulation"""
b2.store('Initialization')

print("Loading I0 for target average firing rate %.2f +/- %.2f Hz" % (Nu_0/b2.Hz,DeltaNu_0/b2.Hz))
Iext_solution = LoadMatrixI0(J,sigma,I0Folder,SimuIdentifier="MatrixI0JSigma" + I0Folder[3:-1] + str(int(ParameterList[5])) + "Hz"+ "_N" + str(N))*IextUnits

t3 = time.time()
print("\tFull Setup took %.3f sec" % (t3-t1))
print("################################################################################")

"""
Run full simulation
I have the I0 needed to have a fixed population firing rate. Now I clear the RateMonitor and create also all the monitors I needed (After running the network a transient time)
"""

"""Restore intial state of the network. we only gonna change the I0 current""" 
b2.restore('Initialization')

"""Select Iext = Iext_solution""" 
Network.Iext = Iext_solution

"""Print simulation parameters"""
print("I0 = %2.4f, J = %.4f mV, sigma = %.2f mV, N = %d, alpha = %.1f, InitialConditions = %s with %d clusters" % (Iext_solution/IextUnits, J/b2.mV, sigma/b2.mV, N, alpha, mode, clusters))
"""Run network for a transient time"""
print("Running transient of full simulation...")
b2.run(transientTime)

print("Running the first part of the full simulation...\n")
FirstsimulationTime = simulationEnd-transientTime-SavingTime
b2.run(FirstsimulationTime)

"""
Setting monitors

Parameters
----------
RateMon : Population rate monitor
VoltageMon : Monitor of membrane voltage for all neurons. Needed for calculate the synchrony index
SpikeMon : spike monitor for all the population. Needed to make the raster plot
VarMon : Monitoring the full variables for equispaced 20 neurons in the population
"""

print("Setting monitors...")
RateMon = b2.PopulationRateMonitor(Network)
VoltageMon = b2.StateMonitor(Network, 'v', record=True)
SpikeMon = b2.SpikeMonitor(Network, record=True)
VarMon = b2.StateMonitor(Network, variables=True, record = np.arange(0, N, N/Neurons2Save))

print("Running the saving part of the full simulation...\n")
b2.run(SavingTime, report='stdout', report_period=30*b2.second)

Network_mean_rate, RateVsTime = AvgPopulationRate(RateMon,simulationEnd-SavingTime,VectorRateMode=True)
print("Avg. Firing Rate = %.3f Hz" % Network_mean_rate)

"""
Create and formating the output of the simulation

Parameters
----------
ISI_times : Interspike interval times
IndISI_times : Interspike interval times for a single neuron
RasterMatrix : Matrix of spiks vs times, for the raster plot
MeanIndivRate : Mean individual firing rate of the neurons, obtained as the nummber of spikes divided by the number of neurons and the simulation time
RateMatrix : Smoothed population rate vs time
FullMonitorMatrix : Matrix with the full variables for selected 20 neurons

"""

ISI_times, IndISI_times = CalculateISI(SpikeMon,transientTime/b2.ms)
RasterMatrix = np.transpose([SpikeMon.t/b2.ms, SpikeMon.i]) # Save the spiking neuron index i, and the spiking time of that neuron
MeanIndivRate = np.mean(SpikeMon.count/SavingTime)/b2.Hz
RateMatrix = np.transpose([VarMon.t/b2.ms, RateVsTime])


if ChosenModel==0 or ChosenModel==1 or ChosenModel==5:
	FullMonitorMatrix = np.transpose(np.concatenate(([VarMon.t/b2.ms], VarMon.v/b2.mV, VarMon.u/uUnits, VarMon.RecurrentInput/b2.mV)))
	NumberOfColumns = 4
elif ChosenModel==2:
	FullMonitorMatrix = np.transpose(np.concatenate(([VarMon.t/b2.ms], VarMon.v/b2.mV, VarMon.RecurrentInput/b2.mV)))	
	NumberOfColumns = 3
elif ChosenModel==3:
	FullMonitorMatrix = np.transpose(np.concatenate(([VarMon.t/b2.ms], VarMon.v/b2.mV, VarMon.w/wUnits, VarMon.RecurrentInput/b2.mV)))
	NumberOfColumns = 4
elif ChosenModel==4:
	FullMonitorMatrix = np.transpose(np.concatenate(([VarMon.t/b2.ms], VarMon.v/b2.mV, VarMon.RecurrentInput/b2.nA)))
	NumberOfColumns = 3
SynchIndex, AverageMembVoltage = CalculSynchInd(VoltageMon,fullResult=True) # Calculate Synchrony Index
BinderCumul = CalculBinderCumul(VoltageMon) # Calculate Binder Cumulant 

AverageMembMatrix = np.transpose([VarMon.t/b2.ms, AverageMembVoltage])

""" 
Generate a summary of the simulation
Note that the Mean Firing Rate of the neurons is gonna be of the order of the average firing rate of the population, the stochastic population oscillator frequency should be obtained by doing the autocorrelation of the population firing rate across time
"""
SimulationSummary = "Simulation Summary\n I0 = {I0} mV\n J = {J} mV\n sigma = {sigma} mV\n N = {N}\n alpha = {alpha}\n InitialConditions = {mode} with {cluster} clusters\n Synchrony Index = {SynchIndex}\n Binder Cumulant = {BinderCumul}\n Connection Probability = {epsilon}\n Avg. Population Firing Rate = {AvRate} Hz\n Avg. Individual Neuron Firing Rate = {MeanIndivRate} Hz"\
					.format(I0=Iext_solution/b2.mV, J=J/b2.mV, sigma=sigma/b2.mV, N=N, alpha=alpha, mode=mode, cluster=clusters, SynchIndex=SynchIndex, BinderCumul = BinderCumul, epsilon=epsilon, AvRate=Network_mean_rate, MeanIndivRate=MeanIndivRate)

"""Create a string which identifies the simulation, in order to be able to save the .csv files and the figures"""
SimuIdentifier = "J=" + "{:.1f}".format(J/b2.mV) + "_sigma=" + "{:.4f}".format(sigma/b2.mV) + "_N=" + str(N) + "_alpha=" + str(alpha)+ "_nu0=" + str(Nu_0/b2.Hz) + '_InitCond=' + str(mode) + str(clusters)#+"_"+ExcInhTheta

print("Saving text files...")
"""Save ISI"""
ISIFile = WorkingDirectory + CSVFolder + "ISI_" + SimuIdentifier +".csv"
header="Interspike Intervals"
np.savetxt(ISIFile,ISI_times, delimiter = " ",header=header,fmt="%.6f")
"""Save Raster"""
RasterFile = WorkingDirectory + CSVFolder + "Raster_" + SimuIdentifier +".csv"
header = "t \t Spike time"
np.savetxt(RasterFile,RasterMatrix, delimiter = " ",header=header,fmt="%.6f")
"""Save Population Rate vs time"""
RateFile = WorkingDirectory + CSVFolder + "PopRate_" + SimuIdentifier +".csv"
header = "t \t Population Rate"
np.savetxt(RateFile,RateMatrix, delimiter = " ",header=header,fmt="%.6f")
"""Save Average Membrane Voltage vs time"""
AverageMembFile = WorkingDirectory + CSVFolder + "AvMembVolt_" + SimuIdentifier +".csv"
header = "t \t Average Membrane Voltage"
np.savetxt(AverageMembFile,AverageMembMatrix, delimiter = " ", header=header,fmt="%4.2f %2.6f")
"""Save Simulation Summary"""
SummaryFile = WorkingDirectory + CSVFolder + "Summary_" + SimuIdentifier +".csv"
with open(SummaryFile, "w") as text_file:
    text_file.write(SimulationSummary)
"""Save full variables for Neurons2Save selected neurons"""
FullMonitoringFile = WorkingDirectory + CSVFolder + "FullMonitor_" + SimuIdentifier +".csv"
header = "".join([" t "]+ [" v_%d " % i for i in range(0,Neurons2Save)] + [" u_%d " % i for i in range(0,Neurons2Save)]  + [" Inotnoisy_%d " % i for i in range(0,Neurons2Save)]) 
np.savetxt(FullMonitoringFile,FullMonitorMatrix, delimiter = " ", header = header,fmt="%4.2f"+ " %2.4f"*Neurons2Save*(NumberOfColumns-1))

t4 = time.time()
print("\tFull Simulation took %.3f sec" % (t4-t1))
print(SimulationSummary)

"""Create and save plot"""
#PlotFile = WorkingDirectory + FiguresFolder + SimuIdentifier +".png"
#CreatePlot2DModel(simulationEnd,ISI_times,SpikeMon,VarMon,Iext_solution/b2.mV,PlotFile,IntervalLong=PlotingTime,NeuronIndex=0,v_th=v_th,bin_lim = None, num_bins = 15)
PlotFile2 = WorkingDirectory + FiguresFolder + SimuIdentifier + ".svg"
# Izhi Type I
#CreatePlot2DModelSkinner(simulationEnd,ISI_times,SpikeMon,VarMon,Iext_solution/IextUnits,PlotFile2,IntervalLong=PlotingTime,NeuronIndex=0,v_th=v_th,bin_lim = None, num_bins = 15)
# Izhi Type II
binLim=None
CreatePlot2DModel(simulationEnd,ISI_times,SpikeMon,VarMon,PlotFile2,N,IntervalLong=PlotingTime,NeuronIndex=0,v_th=v_th,bin_lim = binLim, num_bins = 30, Theta=ThetaInput, ChosenModel=ChosenModel)
print('Done!')
print("################################################################################")
