import brian2 as b2
import numpy as np
from SimulationParametersVarDrive import *
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
Define the new biexponential synapse model. Same for all simulated networks 

Parameters
----------
s : postsynaptic current
x : auxiliary variable
DeltaX : increment of the variable x after a spike
tau_d : decay time of the postsynaptic current
tau_r : rise time of the postsynaptic current

"""
tau_fall = 2.*b2.ms
tau_rise = 0.3*b2.ms
gms = 1.65*b2.nS
delays = 0.8*b2.ms
# Normalizing factor for the IPSC conductance height
c_fall = 1./tau_fall; c_rise = 1./tau_rise

norm_syn = 1./( np.exp(-c_fall*np.log(c_rise/c_fall)/(c_rise-c_fall)) - np.exp(-c_rise*np.log(c_rise/c_fall)/(c_rise-c_fall)) )


Tot_eqsVia = '''
                     Isyn = -(v-Esyn)*norm_syn*(g_fall-g_rise)/N : ampere
                     dg_fall/dt = -g_fall/tau_fall : siemens
                     dg_rise/dt = -g_rise/tau_rise : siemens
                     Esyn : volt
             ''' 
synaps_eqsVia = "DeltaX : siemens"
synaps_actionVia = 'g_fall_post+=DeltaX;g_rise_post+=DeltaX' #

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

ExcInput = "I0 = Iext : volt"+"\n""Iext : volt"
#ExcInput = "I0 = Iext*(1.+1.*sin(2*pi*f*t)) : volt"+"\n""Iext : volt"

if ConductanceBased:
	TotalInput = " + I0/tau_m + sigma/sqrt(tau_m)*xi : volt" + "\n" + Tot_eqsVia +"\n"+"stot : 1"+"\n"+ExcInput     
else:
	TotalInput = " + I0/tau_m + RecurrentInput/tau_m + sigma/sqrt(tau_m)*xi : volt" + "\n" + "RecurrentInput = - J/N * stot : volt"+"\n"+"stot : 1"+"\n"+ExcInput

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

IzhiTypeII = "dv/dt = 0.04/mV/ms * v*v + 5/ms *v + 140*mV/ms - u/pF" + TotalInput + "\n" + "du/dt = aIzTyII*(bIzTyII*v-u) : ampere"

#IzhiTypeII = "dv/dt = 0.04/mV/ms*1000 * v*v + 5/ms*1000 *v + 140*mV/ms*1000 - u/uF" + TotalInput + "\n" + "du/dt = aIzTyII*(bIzTyII*v-u) : ampere"

threshold_condIzhiTypeII = "v > v_thIzTyII"
reset_condIzhiTypeII = "v = v_rIzTyII; u += dIzTyII"

v0IzTyII, u0IzTyII, v_thIzTyII, v_rIzTyII = -60.*b2.mV, -15*b2.pA, 30.*b2.mV, -60.*b2.mV

aIzTyII, bIzTyII, cIzTyII, dIzTyII = 0.1/b2.ms, 0.26*b2.nS, -60.*b2.mV, 0*b2.pA

aIzTyIIMod = 0.1/b2.ms
v0IzTyIIMod, v_rIzTyIIMod = -60*b2.mV, -60*b2.mV,
u0IzTyIIMod = -15.*b2.pA

#CapIzTyII = 1.*b2.pF
CapIzTyII = 1.*b2.uF

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

v_rAdExAH = -58.*b2.mV
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
                     betam = km2*exp(v/sigm2)/ms : Hz
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

# Finally the capacitive and clamp currents

if ConductanceBased:
	currs = " ( INa + IKv3 + IKv1 + IL + Isyn)/Cap "
else:
	currs = " ( INa + IKv3 + IKv1 + IL )/Cap "

volt_eqs = "\n dv/dt = " + currs + TotalInput

base_eqs = gat_vars + param_eqs

GuillemModel = base_eqs + volt_eqs

threshold_condGuillem = "v > v_th"
eps = 1e-9*b2.mV

##########################################################################
### LOADING PARAMETERS AND SETTING UNITS
# Set fixed and homogeneous neuron parameters
ENa=+50.*b2.mV; EK=-90.*b2.mV;
km2=.1; kh1=.012; kh2=.2; kn2=.001; ka2=.02;  kn1 = 1.;
sigm1=4.*b2.mV; sigh2=3.5*b2.mV; sign1=12.*b2.mV; sigm2=-13.*b2.mV;
sigh1=-20.*b2.mV; sign2=-8.5*b2.mV; siga1=12.*b2.mV; siga2=-80.*b2.mV;

# Tuned parameters with heterogeneity and properties as observed in electrophysiological recordings
# Load neuron parameters that can be heterogeneous (need to set folder and file names)
params_act = np.transpose(np.loadtxt( "params/intrinsic/active/params_actRoman.dat" )) # Heterogeneous active parameters: Voltage-gated current peak conductances and thetas
folder = "original/";
gLs = np.loadtxt( "params/intrinsic/passive/gLs/" + folder + "gLsRoman.dat" ) # Leakage conductances
ELs = np.loadtxt( "params/intrinsic/passive/ELs/" + folder + "ELsRoman.dat" ) # Leakage reversal potentials
tauvs = np.loadtxt( "params/intrinsic/passive/taums/taumsRomanUniform.dat" ) # Membrane time constants

NeuronModel = ChosenNeuronModelGuillem

Rins = 1.0e+3/gLs
Caps = tauvs/Rins
gNas, gKv3s, gKv1s, thm1s, thh2s, thn1s, tha1s = params_act

# SET HOMOGENEOUS INTRINSIC NEURONAL PARAMETERS IF DESIRED (ONLY FOR ACTIVE OR BOTH FOR ACTIVE AND PASSIVE)
# If considering a homogeneous network, take passive parameters from a single model neuron
gLs, ELs, tauvs, Caps = gLs[NeuronModel], ELs[NeuronModel], tauvs[NeuronModel], Caps[NeuronModel]
# If considering a network with homogeneous intrinsic neuronal active parameters (from kneu set in the function call)
gNas, gKv3s, gKv1s  = gNas[NeuronModel], gKv3s[NeuronModel], gKv1s[NeuronModel]
thm1s, thh2s, thn1s, tha1s = thm1s[NeuronModel], thh2s[NeuronModel], thn1s[NeuronModel], tha1s[NeuronModel]

# Set units
gNas = gNas*b2.nS; gKv3s = gKv3s*b2.nS; gKv1s = gKv1s*b2.nS; # Peak conductances for voltage-gated currents
thm1s = thm1s*b2.mV; thh2s = thh2s*b2.mV; thn1s = thn1s*b2.mV; tha1s = tha1s*b2.mV;
tauvs = tauvs*b2.ms; gLs = gLs*b2.nS; Caps = Caps*b2.nF; ELs = ELs*b2.mV; # Passive properties

print(gNas,gKv3s,gKv1s,tauvs,gLs,Caps,ELs);
'''
Define Wang-Buzsaki model 1996
'''
v_th=-30.*b2.mV; # Voltage threshold to fire a spike
Cm = 1.*b2.uF # /cm**2
gLWB = 0.1*b2.msiemens
ELWB = -65.*b2.mV
ENa = 55.*b2.mV
EK = -90.*b2.mV
gNaWB = 35.*b2.msiemens
gK = 9.*b2.msiemens

WangBuzsakiModelAux = '''
m = alpha_m/(alpha_m+beta_m) : 1
alpha_m = -((v+35.*mV-eps)/10./mV)/(exp(-(v+35.*mV-eps)/10./mV)-1.)/ms : Hz
beta_m = 4.*exp(-(v+60.*mV)/(18.*mV))/ms : Hz

dh/dt = alpha_h*(1-h)-beta_h*h : 1
alpha_h = 0.21*exp(-(v+58.*mV)/(20.*mV))/ms : Hz
beta_h = 3./(exp(-0.1/mV*(v+28.*mV))+1)/ms : Hz

dn/dt = alpha_n*(1-n)-beta_n*n : 1
alpha_n = -0.03*(v+34.*mV)/(exp(-(v+34.*mV)/(10*mV))-1.)/ms/mV : Hz
beta_n = 0.375*exp(-(v+44.*mV)/(80.*mV))/ms : Hz
'''
Area=9.6e-5

WangBuzsakiModel = "dv/dt = (-gNaWB*m**3*h*(v-ENa)-gK*n**4*(v-EK)-gLWB*(v-ELWB))/Cm+Isyn/Caps" + TotalInput + "\n" + WangBuzsakiModelAux
threshold_condWangBuzsaki = "v > v_th"
