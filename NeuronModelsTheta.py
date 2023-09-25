"""
Define the biexponential synapse model. Same for all simulated networks 

Parameters
----------
s : postsynaptic current
x : auxiliary variable
DeltaX : increment of the variable x after a spike
tau_d : decay time of the postsynaptic current
tau_r : rise time of the postsynaptic current

"""
synaps_eqs = '''
				 ds/dt = (-s + x)/tau_d : 1 (clock-driven)
				 DeltaX : 1
				 stot_post = s : 1 (summed) 
				 dx/dt = -x/tau_r : 1 (clock-driven)
             ''' 
synaps_action = 'x += DeltaX' # This statement says after a spike, increase x in postsynaptic neuron by DeltaX 

"""
Define the total input to a neuron, comprising the external "current" (measured in mV as in Brunel 2006), the total recurrent input and the noisy current.

Parameters
----------

tau_m : membrane time constant, in ms. Actually this is a scale factor to compare models with the LIF one. Its values is fixed at 10 ms.
Iext : external DC input, in mV. 
For Theta-Nested Iext*(1.+1.*sin(2*pi*f*t)) Test2: Iext*(0.5+0.5*sin(2*pi*f*t))
J : total synaptic strength/coupling strength, in mV
N : number of neurons in the network
stot : total postsynaptic current, without units. The synapses are current-based.
sigma : noise amplitude, in mV

"""
ExcInhTheta = "Exc" # Depending on wheter Theta input is above or below background input

factorI0 = 1./1.

ExcInput = "I0 = 0.*mV*int(t<=1500*ms)+factorI0*Iext/2.*(1.+sin(2*pi*f*(t-1500*ms)-pi/2))*int(t>1500*ms)*int(t<=2500*ms)+0.*mV*int(t>2500*ms) : volt"+"\n""Iext : volt"
InhInput = "I0 = Iext*int(t<=1500*ms)+Iext/2.*(1.+sin(2*pi*f*(t-1500*ms)+pi/2))*int(t>1500*ms)*int(t<=2500*ms)+Iext*int(t>2500*ms) : volt"+"\n""Iext : volt"
IntermediateInput = "I0 = 0.*mV*int(t<=1500*ms)+Iext*int(t>1500*ms)*int(t<=2500*ms)+0.*mV*int(t>2500*ms) : volt"+"\n""Iext : volt"


TotalInput = " + factorI0*I0/tau_m + RecurrentInput/tau_m + sigma/sqrt(tau_m)*xi : volt" + "\n" + "RecurrentInput = - J/N * stot : volt"+"\n"+"stot : 1"+"\n"+ExcInput
	     

"""
Define the equations of the Izhikevich resonator Type II model. Some parameters are settled in SimulationParameters.py (sigma, J, N)

Parameters
----------

v : membrane voltage, in mV
u :  adaptation variable, in mV
xi : white gaussian noise, in sqrt(s)
v_rIzTyII : after-spike reset voltage, in mV
v_thIzTyII : voltage cut-off value, in mV
aIzTyII : recovery variable time scale, in 1/ms
bIzTyII : sensitivity of u to subthreshold fluctuations of v, no units
cIzTyII : voltage reset, in mV
dIzTyII : recovery variable after-spike adaptation, in mV

"""

IzhiTypeII = "dv/dt = 0.04*nS/mV/pF * v*v + 5*nS/pF *v + 140*mV*nS/pF - u/pF" + TotalInput + "\n" + "du/dt = aIzTyII*(bIzTyII*v-u) : ampere"
#-1/5*5 1/23*140

threshold_condIzhiTypeII = "v > v_thIzTyII"
reset_condIzhiTypeII = "v = v_rIzTyII; u += dIzTyII"

v0IzTyII, u0IzTyII, v_thIzTyII, v_rIzTyII = -60.*b2.mV, -15*b2.pA, 30.*b2.mV, -60.*b2.mV

aIzTyIINorm, bIzTyII, cIzTyII, dIzTyII = 0.1/b2.ms, 0.26*b2.nS, -60.*b2.mV, 0*b2.pA

aIzTyIIMod = 0.1/b2.ms
v0IzTyIIMod, v_rIzTyIIMod = -20*b2.mV, -20*b2.mV,
u0IzTyIIMod = -15.*b2.pA

CapIzTyII = 1.*b2.pF

uUnitsIzTyII = 1.*b2.nA


"""
Define the equations of the Izhikevich FS PV+ Skinner Type I model. Some parameters are settled in SimulationParameters.py (sigma, J, N)

Parameters
----------

v : membrane voltage, in mV
u :  adaptation variable, in pA
v_rIzTyI : after-spike reset voltage, in mV
v_thIzTyI : voltage cut-off value, in mV
CmIzTyI : membrane capacitance, in pF
vtIzTyI : instantaneous threshold potential, in mV
aIzTyI : recovery variable time scale, in 1/ms
bIzTyI : sensitivity of u to subthreshold fluctuations of v, in nS
cIzTyI : voltage reset, in mV
dIzTyI : recovery variable after-spike adaptation, in pA
k : scaling factor used to adjust the spike width after the instantaneous threshold vtIzTyI

"""

IzhiTypeI = "k = k_low*int(v <= vtIzTyI) + k_high*int(v > vtIzTyI) : siemens/volt" + "\n" + "dv/dt = (k*(v-v_rIzTyI)*(v-vtIzTyI) - u )/Cm" + TotalInput +"\n" + "du/dt = aIzTyI*(bIzTyI*(v-v_rIzTyI)-u) : ampere"
		   
threshold_condIzhiTypeI = "v > v_thIzTyI"
reset_condIzhiTypeI = "v = v_rIzTyI; u += dIzTyI"

v0IzTyI, u0IzTyI, v_thIzTyI, v_rIzTyI, vtIzTyI, CmIzTyI, k_low, k_high = -60.*b2.mV, -0.5*b2.nA, 2.5*b2.mV, -60.6*b2.mV, -43.1*b2.mV, 90*b2.pF, 1.7*b2.nS/b2.mV, 14.*b2.nS/b2.mV

aIzTyI, bIzTyI, cIzTyI, dIzTyI = 0.1/b2.ms, -0.1*b2.nS, -67*b2.mV, 0.1*b2.pA

uUnitsIzTyI = 1.*b2.pA

"""
Define the equations of the LIF model. The parameters are settled in SimulationParameters.py

Parameters
----------

v : membrane voltage, in mV
tau_m : membrane time constant, in ms
Iext : external DC input, in mV. For Theta Drive I0*(0.7+0.3*sin(2*pi*f*t))
J : Total synaptic strength/coupling strength, in mV
N : number of neurons in the network
stot : total postsynaptic current, without units. The synapses are current-based.
sigma : noise amplitude, in mV
xi : white gaussian noise, in sqrt(ms)
v_r : after-spike reset voltage, in mV
v_th : voltage threshold, in mV

"""
LIFModel = "dv/dt = -v/tau_m" + TotalInput

threshold_condLIF = "v > v_thLIF"
reset_condLIF = "v = v_rLIF"

v0LIF, v_thLIF, v_rLIF = 16.*b2.mV, 20.*b2.mV, 14.*b2.mV

"""
Define the equations of the AdEx model. The parameters are settled in SimulationParameters.py

Parameters
----------

v : membrane voltage, in mV
tau_m : membrane time constant, in ms
Iext : external DC input, in mV. For Theta-Nested I0*(0.7+0.3*sin(2*pi*f*t)) Test2: I0*(0.5+0.5*sin(2*pi*f*t))
J : Total synaptic strength/coupling strength, in mV
N : number of neurons in the network
stot : total postsynaptic current, without units. The synapses are current-based.
sigma : noise amplitude, in mV
xi : white gaussian noise, in sqrt(s)
v_r : after-spike reset voltage, in mV
v_th : voltage threshold, in mV

"""

AdExModel = "dv/dt = -glAdEx/CAdEx*(v-ElAdEx) + glAdEx*DeltaTAdEx/CAdEx*exp((v-V_TAdEx)/DeltaTAdEx)-w/CAdEx" + TotalInput +"\n" + "dw/dt = aAdEx/tau_wAdEx*(v-ElAdEx)-w/tau_wAdEx : ampere"

threshold_condAdEx = "v > v_thAdEx"
reset_condAdEx = "v = v_rAdEx; w += bAdEx"

v_thAdEx = 0.*b2.mV
v_rAdExSNTypeI = -50.*b2.mV
v_rAdExSNTypeII = -38.*b2.mV
tau_wAdExSN, aAdExSN, bAdExSN, CAdExSN, glAdExSN, ElAdExSN, V_TAdExSN, DeltaTAdExSN = 16.*b2.ms, 1.8*b2.nS, 0.*b2.pA, 59.*b2.pF, 2.9*b2.nS, -62.*b2.mV, -42.*b2.mV, 3.*b2.mV

v_rAdExAH = -65.*b2.mV
tau_wAdExAH, aAdExAH, bAdExAH, CAdExAH, glAdExAH, ElAdExAH, V_TAdExAH, DeltaTAdExAH = 50.*b2.ms, 8.0*b2.nS, 0.*b2.pA, 150.*b2.pF, 10.*b2.nS, -58.*b2.mV, -47.5*b2.mV, 0.5*b2.mV 

v0AdEx, w0AdExSN, w0AdExAH = -60.*b2.mV, .01*b2.nA, .4*b2.nA

wUnitsAdEx = 1.*b2.nA

"""
Define the equations of the Via's model. The parameters are settled in SimulationParameters.py

Parameters
----------

"""
# Membrane-potential dynamics equations
gat_vars =       """
                     alpham = -((v-thm1-eps)/sigm1)/(exp(-(v-thm1-eps)/sigm1)-1.)/ms : Hz
                     betam = km2*exp(-v/sigm2)/ms : Hz
                     dm/dt = alpham*(1.-m) - betam*m : 1

                     alphah = kh1/exp(-v/sigh1)/ms : Hz
                     betah = -kh2*(v-thh2)/(exp(-(v-thh2)/sigh2)-1.)/ms/mV : Hz
                     dh/dt = alphah*(1.-h) - betah*h : 1

                     alphan = -(v-thn1)/(exp(-(v-thn1)/sign1)-1.)/ms/mV : Hz
                     betan = kn2/exp(-v/sign2)/ms : Hz
                     dn/dt = alphan*(1.-n) - betan*n : 1

                     alphaa = -(v-tha1)/(exp(-(v-tha1)/siga1)-1.)/ms/mV : Hz
                     betaa = ka2/exp(-v/siga2)/ms : Hz
                     da/dt = alphaa*(1.-a) - betaa*a : 1

                     INa = - gNa*m*m*m*h*(v-ENa) : amp
                     IKv3 = - gKv3*n*n*n*n*(v-EK) : amp
                     IKv1 = - gKv1*a*a*a*a*(v-EK) : amp
                     IL = - gL*(v-EL) : amp
                 """

param_eqs =   """
                  I_sin : amp 
                  g_sin : siemens 
                  EL : volt
                  tauv : second
                  gL : siemens
                  Cap : farad
                  gNa : siemens
                  gKv3 : siemens
                  thm1 : volt
                  thh2 : volt
                  thn1 : volt
                  gKv1 : siemens 
                  tha1 : volt 
              """

# External currents
ext_eqs = TotalInput

# Finally the capacitive and clamp currents
currs = " ( INa + IKv3 + IKv1 + IL )/Cap "

volt_eqs = "\n dv/dt = " + currs + TotalInput

base_eqs = gat_vars + param_eqs

GuillemModel = base_eqs + volt_eqs

threshold_condGuillem = "v > v_th"

params = np.loadtxt("params.dat"); # Parameters: EL, Rin and tauv, gNaf, rat (=gKf/gNaf), thm1, thh2, thn1, and gKv1f, tha1

ELs, Rins, tauvs, gNafs, rats, thm1s, thh2s, thn1s, gKv1fs, tha1s = params[:,ChosenNeuronModelGuillem]

# Providing parameters with units
# Effective, membrane-integrated, capacitances can be obtained from time constants and input resistances. One for each neuron in the network.

#ELs = -40
print(ELs)
ELs = ELs*b2.mV
Rins = Rins*b2.Mohm
tauvs = tauvs*b2.ms
thm1s, thh2s, thn1s, tha1s = thm1s*b2.mV, thh2s*b2.mV, thn1s*b2.mV, tha1s*b2.mV

# Setting parameters

v_th=-30.*b2.mV; # Voltage threshold to fire a spike

ENa=+50.*b2.mV; EK=-90.*b2.mV; eps = 1.e-7*b2.mV;

km2=.1; kh1=.012; kh2=.2; kn2=.001; ka2=.02;

sigm1=4.*b2.mV; sigh2=3.5*b2.mV; sign1=12.*b2.mV; sigm2=-13.*b2.mV; sigh1=-20.*b2.mV; sign2=-8.5*b2.mV; siga1=12.*b2.mV; siga2=-80.*b2.mV;

# Capacitance
Caps = tauvs/Rins; gLs = 1./Rins;

# Active conductances
# Heterogeneous peak conductances for active currents. Above are loaded normalized to Golomb values times capacitance.
gNas=np.array(gNafs)*112.5e+3*b2.Hz*Caps;
gKv3s=np.array(gNafs)*np.array(rats)*225.e+3*b2.Hz*Caps;
gKv1s=np.array(gKv1fs)*225.e+3*b2.Hz*Caps;
