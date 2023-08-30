"""File containing all the auxiliar functions needed for the simulation

"""
import numpy as np
import brian2 as b2
import os
import sys
import glob
import matplotlib.pyplot as plt
import scipy
from NeuronModelsVarDrive import *
from SetupModels import *

def LoadMatrixI0(J,sigma,I0Folder,SimuIdentifier="/MatrixI0JSigma"):
	"""
	Function for loading the constant input I0 needed to maintain the average firing rate constant

	Parameters
	----------
	J : Coupling strength
	sigma : Noise intensity
	N : population size

	"""
	WorkingDirectory = os.getcwd() + I0Folder
	Matrix = np.loadtxt(WorkingDirectory + SimuIdentifier)
	IndexJ = np.argwhere(abs(Matrix[:,1]-J/b2.mV)<1e-5)
	IndexSigma = np.argwhere(abs(Matrix[:,2]-sigma/b2.mV)<1e-4)
	Indexes = np.intersect1d(IndexJ,IndexSigma) #Intersection between two arrays
	if Indexes.size == 0:
		raise ValueError	
	I0 = Matrix[Indexes,0]
	return I0
##################################################################################
def SaveI0SigmaJ(I0Folder,alpha,nu,mode,clusters):
	"""
	Function to load all the I0 files, and saving all the info as a matrix in a single file

	Parameters
	----------
	I0Folder : Name of the folder where the files eith the I0 information are saved
	alpha : Synaptic time scale
	mode : Initialization mode: "cluster" or "random"
	clusters : number of initialization clusters
	nu : Desired Average firing rate

	Returns
	-------
	None, but saves in I0folder a file called MatrixI0JSigma with all the information contained in the I0_*.txt files

	Additional
	----------	
	Not implemented as a function because the parallel arrays would imply that the same file is gonna be written over an over again, e very time the function is called
	"""
	WorkingDirectory = os.getcwd() + I0Folder

	alpha = 1
	mode = "cluster"
	clusters = 1
	nu = 30

	SimuIdentifier = "*_alpha=" + str(alpha) + "*_nu0=" + str(nu) + '*_InitCond=' + str(mode) + str(clusters) + "*.csv"

	## Find and open the files
	os.chdir(WorkingDirectory)
	ListsOfLists = []

	NumberOfSimulations = len(glob.glob(SimuIdentifier))
	MatrixData = []

	for file1 in glob.glob(SimuIdentifier):
		A = np.loadtxt(file1)
		A[1] = A[1]*1000
		A[2] = A[2]*1000
		MatrixData.append(A)
	MatrixData = np.row_stack(MatrixData)	
	np.savetxt("MatrixI0JSigma",MatrixData)
##################################################################################
def SetInitialVoltages(N, mode='random', clusters=1, v_th=30., v_r=-60., v0=-60., u0 =-15.):
	"""
	Function for setting the initial conditions. Initialization could be random or cluster (equally spaced clusters)

	Parameters
	----------
	N : number of neurons
	mode : Initialization, could be 'random','cluster'
	clusters : number of clusters
	v_th : voltage threshold
	v_r : voltage reset
	v0 : initial voltage for 1 cluster
	u0 : initial recovery variable for 1 cluster

	Returns
	-------

	V : N-length array of initial membrane potential
	U : N-length array of initial recovery variable. For 1-D models as LIF, you can use this output as a dummy variable

	"""
	if mode.lower() == 'random':
		V = (v_th-v_r)/b2.mV*np.random.uniform(size=N)+v_r/b2.mV
		U = np.random.uniform(u0+3.*b2.pA, u0-3.*b2.pA, size=N)
	elif mode.lower()=='cluster':
		DeltaV = (v_th-v_r)/clusters
		DeltaU = 1.
		V = np.concatenate([([v0+(-1)**i*DeltaV*i]*int(N*1./clusters)) for i in range(clusters)], axis=0)
		U = np.concatenate([([u0+(-1)**i*DeltaU*i]*int(N*1./clusters)) for i in range(clusters)], axis=0)
	if np.shape(V)[0]==N:	return V, U
	if np.shape(V)[0]<N:	return np.append(V,[v0]*(N-clusters*int(N*1./clusters))), np.append(U,[u0]*(N-clusters*int(N*1./clusters)))
##################################################################################
def RunLite(I0_init, verbose=False, transientTime=1200*b2.ms, testTime=800*b2.ms):
	"""
	Function to run a lighter simulation, only to find the mean firing rate and adjust the input I0

	Parameters
	----------
	I0_init : DC input
	verbose : True for printing out the simulation report. False (default)
	transientTime : Time needed for the network to reach a stationary regime, in ms. This transient is not of our interest for now.
	testTime : time used to calculate the population average frequency. Around the transition the fluctuation could be high, so we need a larger amount of time to determine the average frequency in a feasible and repetitive way
	
	Returns
	-------
	Nu : Average firing rate, in Hz

	"""
	b2.restore('Initialization') # Restore initial network state
	
	Network.Iext = I0_init

	print("Running transient...")
	b2.run(transientTime) # Run network for a transient time
	## Run the network for a testTime to calculate the mean firing rate
	if verbose:
		b2.run(testTime, report='stdout', report_period=30*b2.second)
	else:
		print("Running simulation in lite mode...")
		b2.run(testTime)
		
	Network_mean_rate = AvgPopulationRate(RateMon,transientTime) # Calculate the mean firing rate
	print("Avg. Firing Rate = %.3f spk/s\n" % Network_mean_rate)
	Nu = Network_mean_rate*b2.Hz
	return Nu
##################################################################################
def BrunelSelfConsistentCurrent(J=10., sigma=1., nu=30., Vt=20., tau_m=10.):
	"""
	Using the approximate relation between the network parameters and the average frequency we obtain the DC current necessary to maintain the desired frequency.
	It works as an initial guess for the I0 current for the simulations. It should be adjusted running the bisectionSearch function.
	Approximation strictly valid for (Vt+J*nu*taum-I0)/sigma >> 1, but it works pretty good for other regions.

	Parameters
	----------
	J : coupling strength
	sigma : noise amplitude
	nu : desired average frequency
	Vt : voltage threshold
	tau_m : membrane time constant

	Returns
	-------
	rootGuess : Approximated I0 current

	"""
	nu = nu/1000. # To convert nu to 1/ms
	yt = lambda I0 : (Vt+J*nu*tau_m-I0)/sigma
	
	f = lambda I0 : yt(I0)/np.sqrt(np.pi)*np.exp(-yt(I0)**2) # Approximation of the integral
	f2 = lambda I0 : abs(nu*tau_m - f(I0))

	Iseed = np.linspace(10,60,100)
	index = 0
	xopt = np.zeros(np.shape(Iseed))

	for seed in Iseed:
		xopt[index] = fmin(f2,np.array([seed]),disp=0)
		index += 1

	fEval = f2(xopt)	
	minimum = np.argmin(fEval)
	rootGuess = xopt[minimum] 
	return rootGuess # Approximated root
##################################################################################
def bisectionSearch(I0_init=0.*b2.mV,I0_end=1.*b2.mV,Nu_0=30.*b2.Hz, DeltaNu_0=0.5*b2.Hz,maxIter=25,I0_lowbound=-10.*b2.mV,I0_uppbound=10.*b2.mV,I0Units=1.*b2.mV,transientTime=1200*b2.ms, testTime=800*b2.ms):
	"""
	Function to do a bisection search to adjust I0 in order to obtain the average firing rate between Nu_0-DeltaNu_0 and Nu_0+DeltaNu_0

	Parameters
	----------
	I0_init : Initial parameter to find the solution
	I0_end : Final parameter to find the solution
	Nu_0 : Desired Frequency
	DeltaNu_0 : Absolute error in the determination of the average frequency
	maxIter : Maximum number of iterations in the bisection algorithm
	I0_lowbound : I0 lower bound. In case the solution is not in the initial interval [I0_init,I0_end], the search restarts but in the interval [I0_lowbound,I0_uppbound]
	I0_uppbound : I0  upper bound.
	I0Units : Units of the input current

	Returns
	-------
	I0_mid : Current which generates average firing frequency around the desired frequency

	"""
	steps = 0
	Network_mean_rate_1 = RunLite(I0_init,transientTime=transientTime, testTime=testTime) # Calculate average firing rate in the extrems of the search interval
	if (abs(Network_mean_rate_1-Nu_0) < DeltaNu_0):
		print("Bisection search converged to %2.12f after %d steps" % (I0_init/I0Units, steps)) 
		return I0_init
	Network_mean_rate_2 = RunLite(I0_end,transientTime=transientTime, testTime=testTime)
	if (abs(Network_mean_rate_2-Nu_0) < DeltaNu_0):
		print("Bisection search converged to %2.12f after %d steps" % (I0_end/I0Units, steps)) 
		return I0_end

	steps+=1	
	I0_mid = (I0_init+I0_end)/2 # Select middle point
	Network_mean_rate_mid = RunLite(I0_mid,transientTime=transientTime, testTime=testTime) # Calculate average firing rate
	## Check whether the middle point is solution
	if (abs(Network_mean_rate_mid-Nu_0) < DeltaNu_0):
		print("Bisection search converged to %2.12f after %d steps" % (I0_mid/I0Units, steps)) 
		return I0_mid
	I0_1 = I0_init
	I0_2 = I0_end
	if (Network_mean_rate_1-Nu_0)/b2.Hz*(Network_mean_rate_2-Nu_0)/b2.Hz>0: 
		print("The ends of the searching interval had the same sign, none root possible. Changing to default extremes")
		I0_1 = I0_lowbound
		I0_2 = I0_uppbound
		Network_mean_rate_1 = RunLite(I0_1,transientTime=transientTime, testTime=testTime)
		Network_mean_rate_2 = RunLite(I0_2,transientTime=transientTime, testTime=testTime)

	while (abs(Network_mean_rate_mid-Nu_0) >= DeltaNu_0 and steps <= maxIter):
		print("Entering loop...")
		steps+=1
		# Decide the side to repeat the steps
		if ((Network_mean_rate_mid-Nu_0)*(Network_mean_rate_1-Nu_0) > 0):
			I0_1 = I0_mid
			Network_mean_rate_1 = RunLite(I0_1,transientTime=transientTime, testTime=testTime)	
		else:
			I0_2 = I0_mid
		# Find middle point
		I0_mid = (I0_1+I0_2)/2
		Network_mean_rate_mid = RunLite(I0_mid,transientTime=transientTime, testTime=testTime)
	if steps > maxIter: raise ValueError
	print("Bisection search converged to %2.12f after %d steps" % (I0_mid/I0Units, steps))
	return I0_mid
##################################################################################
def bisectionSearchBrunel(I0_init,DeltaI= 10*b2.mV,Nu_0=30.*b2.Hz, DeltaNu_0=0.5*b2.Hz,maxIter=20):
	## Run simulation
	steps = 0
	Network_mean_rate = RunLite(I0_init)
	# By trial-and-error we choose a limit for the firing rate. If is bigger than 40 Hz, we select I0_init/2 as an initial point. This is due to the fact that the approximation for the self-consistent equation fails sometimes
	if Network_mean_rate > 40.*b2.Hz: 
		I0_init=I0_init/2
		Network_mean_rate = RunLite(I0_init)
	## Check whether the firing rate is inside the interval (Nu_0-DeltaNu_0, Nu_0+DeltaNu_0)
	if (abs(Network_mean_rate-Nu_0) <= DeltaNu_0): 
		print("Bisection search converged to %2.4f mV after %d steps" % (I0_init/b2.mV, steps))
		return I0_init/b2.mV
	
	if Network_mean_rate > Nu_0+DeltaNu_0: 
		I0_1 = I0_init - DeltaI
		I0_2 = I0_init
		Network_mean_rate_1 = RunLite(I0_1)
		Network_mean_rate_2 = Network_mean_rate
		if (abs(Network_mean_rate_1-Nu_0) < DeltaNu_0):
			print("Bisection search converged to %2.4f mV after %d steps" % (I0_1/b2.mV, steps))
			return I0_1/b2.mV
	elif Network_mean_rate < Nu_0-DeltaNu_0:
		I0_1 = I0_init
		I0_2 = I0_init + DeltaI
		Network_mean_rate_1 = Network_mean_rate
		Network_mean_rate_2 = RunLite(I0_2)
		if (abs(Network_mean_rate_2-Nu_0) < DeltaNu_0):
			print("Bisection search converged to %2.4f mV after %d steps" % (I0_2/b2.mV, steps))
			return I0_2/b2.mV
		
	if (Network_mean_rate_1-Nu_0)/b2.Hz*(Network_mean_rate_2-Nu_0)/b2.Hz>0: 
		print("The ends of the searching interval had the same sign, none root possible. Changing to default extremes")
		I0_1 = 0.*b2.mV
		I0_2 = 50.*b2.mV
		Network_mean_rate_1 = RunLite(I0_1)
		Network_mean_rate_2 = RunLite(I0_2)
	steps+=1	
	I0_mid = (I0_1+I0_2)/2
	Network_mean_rate_mid = RunLite(I0_mid)
	## Check whether the middle point is solution
	if (abs(Network_mean_rate_mid-Nu_0) < DeltaNu_0):
		print("Bisection search converged to %2.4f after %d steps" % (I0_mid/b2.mV, steps)) 
		return I0_mid/b2.mV

	while (abs(Network_mean_rate_mid-Nu_0) >= DeltaNu_0 and steps <= maxIter):
		print("Entering loop...")
		steps+=1
		# Decide the side to repeat the steps
		if ((Network_mean_rate_mid-Nu_0)*(Network_mean_rate_1-Nu_0) > 0):
			I0_1 = I0_mid
			Network_mean_rate_1 = RunLite(I0_1)	
		else:
			I0_2 = I0_mid
		# Find middle point
		I0_mid = (I0_1+I0_2)/2
		Network_mean_rate_mid = RunLite(I0_mid)
	if steps > maxIter: raise ValueError
	print("Bisection search converged to %2.4f after %d steps" % (I0_mid, steps))
	return I0_mid/b2.mV
##################################################################################
def AvgPopulationRate(RateMon,transientTime,VectorRateMode=False):
	"""
	Function for calculate the average firing rate for the population. 

	Parameters
	----------
	RateMon : Rate monitor
	transientTime : Duration of the transitory regime, in ms
	VectorRateMode : Boolean to decide whther return the smoothed population rate across time

	Returns
	-------
	Network_mean_rate : Average population rate of the network
	VectorRates : Smoothed population rate in function of time

	"""
	## Processing (to extract population rates and spike-trains) 
	Rates = [RateMon.t/b2.ms,RateMon.smooth_rate(window='flat',width=1.15*b2.ms)/b2.Hz]
	index_testTime = int(np.where(abs(Rates[0]-transientTime/b2.ms) < 0.0001)[0]) # Get index that corresponds to the end of transientTime	
	VectorRates = Rates[1][index_testTime:]
	Network_mean_rate = np.sum(VectorRates)/len(VectorRates) # Starts at transientTime ms and ends at transientTime+testTime ms
	if VectorRateMode: return Network_mean_rate, VectorRates
	return Network_mean_rate	
##################################################################################
def CalculateISI(SpikeMon,transientTime,neuronInd_index=0):
	"""
	Function to extract the Interspike intervals for the whole population, and for an individual neuron

	Parameters
	----------
	SpikeMon : Spike Monitor
	transientTime : Time elapsed since the beginning of the simulation that is discarded
	neuronInd_index : Index of the neuron we want to extract its ISI

	Returns
	-------
	isi_lst : Interspike intervals list for all the neurons in the network
	isi_ind : Interspike intervals for a single neuron

	"""
	isi_lst = []
	isi_ind = []
	spike_dict = SpikeMon.spike_trains()
	## Loop over all the neurons
	for neuron_index in spike_dict.keys():
		spike_times = spike_dict[neuron_index]/b2.ms
		clean_spike_times = [x for x in spike_times if x >= transientTime] # Removes transientTime of simulation
		for dt in np.diff(clean_spike_times):
			isi_lst.append(dt)
	spike_times = spike_dict[neuron_index]/b2.ms
	clean_spike_times = [x for x in spike_times if x >= transientTime] # Removes transientTime of simulation
	for dt in np.diff(clean_spike_times):
		isi_ind.append(dt)
	return isi_lst, isi_ind
##################################################################################
def CalculSynchInd(VoltageMon,fullResult=False):
	"""
	Function to calculate the population synchrony index,a ccording to Brunel & Hansel 2006

	Parameters
	----------
	VoltageMon : Monitor of the membrane voltage for the whole population
	fullResult : Boolean. Indicates whether return the population mean voltage accross time or not

	Returns
	-------
	SynchIndex : Synchrony Index
	PopulationMeanVoltageAcrossTime : Population mean voltage across time. It could serve as a first order proxy to kind of population LFP (if filtered to 'low' frequencies)

	"""
	MembraneVoltageWholePopulation = VoltageMon.get_states(['v'])['v']/b2.mV
	PopulationMeanVoltageAcrossTime = np.nanmean(MembraneVoltageWholePopulation,axis=1)
	PopulationSquaredStdDev = np.nanmean(PopulationMeanVoltageAcrossTime**2) - np.nanmean(PopulationMeanVoltageAcrossTime)**2
	NeuronSquaredStdDev = np.nanmean(MembraneVoltageWholePopulation**2,axis=0) - np.nanmean(MembraneVoltageWholePopulation,axis=0)**2
	MeanNeuronSquaredStdDev = np.nanmean(NeuronSquaredStdDev)
	SynchIndex = np.sqrt(PopulationSquaredStdDev/MeanNeuronSquaredStdDev)
	if fullResult: return SynchIndex, PopulationMeanVoltageAcrossTime
	return SynchIndex
##################################################################################
def CalculBinderCumul(VoltageMon):
	"""
	Function to calculate the Binder Cumulant. In a critical transition (as the one from synchronous to asynchronous state) at the transition point this quantity is independent of the number of neurons. It allow us to find the transition noise without needing high number of neurons in the simulation

	Parameters
	----------
	VoltageMon : Monitor of the membrane voltage for the whole population

	Returns
	-------
	BinderCumul : Binder Cumulant

	"""
	MembraneVoltageWholePopulation = VoltageMon.get_states(['v'])['v']/b2.mV
	PopulationMeanVoltageAcrossTime = np.mean(MembraneVoltageWholePopulation,axis=1)
	PopulationSquaredSecondCumulant = np.mean(PopulationMeanVoltageAcrossTime**2)
	PopulationSquaredFourthCumulant = np.mean(PopulationMeanVoltageAcrossTime**4)
	BinderCumul = 1 - PopulationSquaredFourthCumulant/(3*PopulationSquaredSecondCumulant**2)
	return BinderCumul
##################################################################################
def CreatePlot2DModel(simulationEnd,ISI_times,SpikeMon,VarMon,PlotFile,NumbNeurons, Nu_0=30,IntervalLong=200,NeuronIndex=0,v_th=20*b2.mV,bin_lim = None, num_bins = None, N=1000, Theta=False,ChosenModel=None,FactorIext=1.,ConductanceBased=False,Timestep=1.):
    """
	Function for creating the figure to plot for 2D neuron model
	It will consist in 6 subfigures. First, the ISI histogram for the whole population. 
	Second, population rate vs time. Third, membrane poential of one neuron. 
	Fourth, phase space in a 2D model. Fifth, input current (I0+Isyn) without noise. 
	Sixth, raster plot.

	Parameters
	----------
	simulationEnd : Final time of the simulation
	IntervalLong : How much time before the end of the simulation we want to plot, in ms (default=200)
	ISI_times : Interspike intervals
	SpikeMon : Spike Monitor
	VarMon : State Monitor
	I0 : Constant external input
	PlotFile : Name of the file for saving the figure, with the corresponding extension
	NeuronIndex : Index of the neuron we want to plot in the membrane potential, and the phase space
	v_th : Voltage threshold
	bin_lim : Limit of the range of plotting for the ISI histogram
	num_bins : Number of bins in the histogram

	Returns
	-------
	None (only plots and save the figure)

    """
    print("Creating figure...")
    plt.figure(figsize=(11,6.5))
    xlim_left = simulationEnd/b2.ms-IntervalLong
    xlim_right = simulationEnd/b2.ms
    ## ISI Histogram
    plt.subplot(231)
    ax1 = plt.gca()
    if (bin_lim == None or num_bins == None) and len(ISI_times)>0:
        ## Calculates optimal bin width & number using the Freedman-Diaconis rule
        optimal_bin_width = 2 * scipy.stats.iqr(ISI_times) / len(ISI_times)**(1./3)
        optimal_num_bins = (max(ISI_times) - min(ISI_times))/optimal_bin_width + 0.5
        ## To avoid numerical errors, we put an upper limit to the number of bins
        if optimal_num_bins < 30:
            l = plt.hist(ISI_times, bins=int(optimal_num_bins), histtype='stepfilled', density=True)
        else:
            l = plt.hist(ISI_times, bins=30, histtype='stepfilled', density=True)
    elif len(ISI_times)>0:
        l = plt.hist(ISI_times, bins=num_bins, histtype='stepfilled', range=(0,bin_lim), density=True)
    else:
        l = plt.hist([])
    plt.title('Histogram')
    plt.xlabel('ISI (ms)')
    plt.ylabel('Counts')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    if max(l[0]) is not None:
        1
    else:
        plt.ylim(0, max(l[0]))
    ## Population firing rate across time
    plt.subplot(232)
    Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
    bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)
    l = plt.hist(SpikeMon.t/b2.ms, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, weights=Nu_0**2/(NumbNeurons*1.0)*np.ones(np.shape(SpikeMon.t/b2.ms)) ) #bar stepfilled
    plt.title("Population Rate")
    plt.xlabel('Time (ms)')
    plt.ylabel('Network Rate (Hz)')
    #plt.ylim(0, 200)
    ax4 = plt.gca()
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    if Theta:
        ax4b = ax4.twinx()
        ax4b.plot(VarMon.t/b2.ms,FactorIext*VarMon.Iext[0]/IextUnits,'-c')
        ax4b.set_ylabel('Theta Input (mV)', color='c')
        ax4b.set_ylim(1.*min(VarMon.I0[0]/IextUnits),1.*max(VarMon.I0[0]/IextUnits))
        ax4b.tick_params(axis='y', colors='c')	
    plt.xlim(xlim_left, xlim_right)

    ## Membrane Potential across time
    plt.subplot(233)
    array = np.asarray(VarMon.t/b2.ms)
    idx1 = (np.abs(array - xlim_left)).argmin()
    idx2 = (np.abs(array - xlim_right)).argmin()
    ax1 = plt.gca()
    ax1.plot(VarMon.t/b2.ms, VarMon.v[NeuronIndex].T/b2.mV, 'b')
    ax1.plot(VarMon.t/b2.ms, v_th/b2.mV*np.ones(np.shape(VarMon.t/b2.ms)), '--g')
    ax1.spines['top'].set_visible(False)
    ax2 = ax1.twinx()
    if ChosenModel==0 or ChosenModel==1 or ChosenModel==5:
        ax2.plot(VarMon.t/b2.ms, VarMon.u[NeuronIndex]/uUnits, 'r')
        ax2.set_ylabel('u', color='r')
        ax2.set_ylim(min(VarMon.u[NeuronIndex]/uUnits)-0.5*abs(min(VarMon.u[NeuronIndex]/uUnits)), max(VarMon.u[NeuronIndex]/uUnits)+0.5*abs(min(VarMon.u[NeuronIndex]/uUnits)))
    elif ChosenModel==3:
        ax2.plot(VarMon.t/b2.ms, VarMon.w[NeuronIndex].T/wUnits, 'r')
        ax2.set_ylabel('w', color='r')
        ax2.set_ylim(min(VarMon.w[NeuronIndex].T/wUnits)-0.5*abs(min(VarMon.w[NeuronIndex].T/wUnits)), max(VarMon.w[NeuronIndex].T/wUnits)+0.5*abs(min(VarMon.w[NeuronIndex].T/wUnits)))
    ax2.tick_params(axis='y', colors='red')	
    ax2.spines['top'].set_visible(False)
    plt.title('Single Neuron Trace')
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('V (mV)')
    ax1.set_xlim(xlim_left, xlim_right)
    ax1.set_ylim(min(VarMon.v[NeuronIndex].T[idx1:idx2]/b2.mV)-1, v_th/b2.mV+1)
    ## Empty subplot or Phase space plot in the case of 2D neuron models
    plt.subplot(234)
    ax1 = plt.gca()
    if ChosenModel==0 or ChosenModel==1 or ChosenModel==5:
        ax1.plot(VarMon.v[NeuronIndex].T[idx1:idx2]/b2.mV, VarMon.u[NeuronIndex].T[idx1:idx2]/uUnits, 'b--')
        ax1.set_ylabel('u (nA)')
    elif ChosenModel==3:
        ax1.plot(VarMon.v[NeuronIndex].T[idx1:idx2]/b2.mV, VarMon.w[NeuronIndex].T[idx1:idx2]/wUnits, 'b--')
        ax1.set_ylabel('w (nA)')
    ax1.set_xlabel('v (mV)')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ## Neuron Input curent without noise
    plt.subplot(235)
    ax3 = plt.gca()
    if ConductanceBased:
        factor=10./np.array(Caps/b2.pF)
        ax3.plot(VarMon.t/b2.ms, FactorIext*VarMon.Iext[0]/IextUnits + factor*VarMon.Isyn[NeuronIndex].T/b2.nA, 'g')
        ax3.plot(VarMon.t/b2.ms, FactorIext*VarMon.Iext[0]/IextUnits + factor*np.mean(VarMon.Isyn[NeuronIndex].T/b2.nA) *np.ones(np.shape(VarMon.t/b2.ms)), 'b')
    else:
        ax3.plot(VarMon.t/b2.ms, FactorIext*VarMon.Iext[0]/IextUnits + VarMon.RecurrentInput[NeuronIndex].T/b2.mV, 'g')
        ax3.plot(VarMon.t/b2.ms, FactorIext*VarMon.Iext[0]/IextUnits + np.mean(VarMon.RecurrentInput[NeuronIndex].T/b2.mV) *np.ones(np.shape(VarMon.t/b2.ms)), 'b')

    if ChosenModel==0: # IzhiTypeI
        factor=10/0.09
        ax3.plot(VarMon.t/b2.ms, factor*0.129283*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
    elif ChosenModel==1: # IzhiTypeII
        factor=10
        ax3.plot(VarMon.t/b2.ms, factor*0.1795*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
        ax3.plot(VarMon.t/b2.ms, factor*0.2625*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
        ax3.plot(VarMon.t/b2.ms, factor*0.4225*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
    elif ChosenModel==2: # LIF
        factor=1
        ax3.plot(VarMon.t/b2.ms, factor*20*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
    elif ChosenModel==3:
        factor=10./np.array(CAdEx/b2.pF)
        if ChosenSubmodelAdEx==2: # AdEx AndronovHopf
            ax3.plot(VarMon.t/b2.ms, factor*165*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
            ax3.plot(VarMon.t/b2.ms, factor*184.8613*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
            ax3.plot(VarMon.t/b2.ms, factor*185.2901*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
        else: # AdEx SaddleNode
            ax3.plot(VarMon.t/b2.ms, factor*86.70821005*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
    elif ChosenModel==4: # Via's
        factor=10./np.array(Caps/b2.pF)
    elif ChosenModel==5: # IzhiTypeII Modified
        factor=10/0.09
        ax3.plot(VarMon.t/b2.ms, factor*0.129283*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
        ax3.set_ylim(-30,50)
        #factor=10
        #ax3.plot(VarMon.t/b2.ms, factor*0.1795*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
        #ax3.plot(VarMon.t/b2.ms, factor*0.2625*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')
        #ax3.plot(VarMon.t/b2.ms, factor*0.4225*np.ones(np.shape(VarMon.t/b2.ms)), 'r--')

    plt.title('Synaptic Current')
    ax3.set_ylabel('Net Current (mV)')
    ax3.set_xlabel('Time (ms)')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.set_xlim(xlim_left, xlim_right)
    ## Raster plot
    plt.subplot(236)
    DummyOut, Index1, Index2 = np.unique(SpikeMon.i,return_index=True,return_inverse=True)
    IndexAux = Index2[:21]
    spike_dict = SpikeMon.spike_trains()
    j=1
    #for i in IndexAux:
    #	plt.scatter(spike_dict[i]/b2.ms, np.ones(len(spike_dict[i]))*j, s=15.0, color ='black', marker = '|')
    #	j+=1
    for key in spike_dict.keys():
        if key%int(N/20) == 0:
            plt.scatter(spike_dict[key]/b2.ms, np.ones(len(spike_dict[key]))*key, s=15.0, color ='black', marker = '|')
    plt.xlim(xlim_left, xlim_right)
    plt.xlabel('Time (ms)')
    ax4 = plt.gca()
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['left'].set_visible(False)
    plt.tick_params(labelleft=False,left=False)
    ## For adjusting the spaces between subplots
    plt.tight_layout()
    ## Save Plot
    print("Saving figure...")
    plt.savefig(PlotFile)
######################################################
def PreviousN(N,NAux):
	pos = np.argwhere(abs(NAux-N*np.ones(np.shape(NAux)))<1e-3)
	if pos[0][0]>0:
		return NAux[pos[0][0]-1]
	else:
		print("Error")
		return NAux[pos[0][0]]

