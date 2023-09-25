"""
File to solve the constant DC input for the Izhikevich type 2 resonator network

Code for running and found the I0 input necesary to maintain a determined average firing rate

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
execfile("SimulationParametersCompleteI0.py")
execfile("AuxiliarFunctions.py")
execfile("NeuronModels.py")
execfile("SetupModels.py")

NAux = [800]
NAux2 = [800,3000]#[800,1000,1400,1800,2200,3000]

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
	print(Model+"\n")
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

	if ChosenModel==0 or ChosenModel==1:
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

"""
Set Rate monitor
I need a rate monitor to calculate the population average firing rate, and to find the DC input that could make the population fire at the desired frequency
"""
RateMon = b2.PopulationRateMonitor(Network)
t2 = time.time() # Calculating initialization time
print("\tInitialization took %.3f sec" % (t2-t1))

"""Store initial state for restoring during bisection search and full simulation"""
b2.store('Initialization')
print("J=" + "{:.1f}".format(J/b2.mV) + " sigma=" + "{:.4f}".format(sigma/b2.mV) + " N=" + str(N) + " alpha=" + str(alpha)+ " nu0=" + str(Nu_0/b2.Hz) + ' InitCond=' + str(mode) + str(clusters))

print("Starting adjustement of Iext for target average firing rate %.2f +/- %.2f Hz" % (Nu_0/b2.Hz,DeltaNu_0/b2.Hz))


SimuIdentifier="MatrixI0JSigma" + I0Folder[3:-1] + str(int(ParameterList[5])) + "Hz"+ "_N" + str( N )
WorkingDirectory = os.getcwd() + I0Folder
Matrix = np.loadtxt(WorkingDirectory + SimuIdentifier)
IndexJ = np.argwhere(abs(Matrix[:,1]-J/b2.mV)<1e-5)
IndexSigma = np.argwhere(abs(Matrix[:,2]-sigma/b2.mV)<1e-4)
Indexes = np.intersect1d(IndexJ,IndexSigma) #Intersection between two arrays
if Indexes.size == 0:
	A = Matrix[:,1:3];
	tree = spatial.KDTree(A)
	Iext_Guess = Matrix[tree.query([J/b2.mV,sigma/b2.mV])[1],0]*IextUnits

	I0_solution = bisectionSearch(I0_init=Iext_Guess-0.1*abs(Iext_Guess),I0_end=Iext_Guess+0.1*abs(Iext_Guess),Nu_0=Nu_0, DeltaNu_0=DeltaNu_0,I0_lowbound=16.*IextUnits,I0_uppbound=25.*IextUnits,I0Units=IextUnits,transientTime=transientTime,maxIter=50)/IextUnits
else:	
	print("Already Calculated")
	quit()

t3 = time.time()
print("\tFull Setup took %.3f sec" % (t3-t1))
print("################################################################################")

SimuIdentifier = "J=" + "{:.1f}".format(J/b2.mV) + "_sigma=" + "{:.4f}".format(sigma/b2.mV) + "_N=" + str(N) + "_alpha=" + str(alpha)+ "_nu0=" + str(Nu_0/b2.Hz) + '_InitCond=' + str(mode) + str(clusters)
WorkingDirectory = os.getcwd() 

print("Saving text files...")
"""Save I0"""
I0File = WorkingDirectory + I0Folder + "I0_" + SimuIdentifier +".csv"
I0Text = "{I0}".format(I0=I0_solution) + " " + "{J}".format(J=J/b2.mV) + " " + "{sigma}".format(sigma=sigma/b2.mV)+ " " + str(N)
with open(I0File, "w") as text_file:
    text_file.write(I0Text)
print("Done!")
