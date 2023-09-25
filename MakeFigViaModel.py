import numpy as np
import matplotlib.pyplot as plt
import glob, os
import matplotlib as mlb

mlb.rcParams['font.serif'] = "Arial"
mlb.rcParams['font.size'] = "10"

#Synapse = "Shunt"
Synapse = "Hyper"

Folder = '../CSV_GuillemVarDrive'+ Synapse +'/*'
SimuIdentifierHyper = '*sigma=10.0110*N=3000*I0=65.*'
SimuIdentifierShunt = '*sigma=7.0110*N=3000*I0=70.*'

SimuIdentifier1=SimuIdentifierHyper
edgesbinHyper = np.arange(0,200,1)
edgesbinShunt = np.arange(0,100,1)
edgesbin = edgesbinHyper

if Synapse =="Shunt":
	ISII0_10 = glob.glob(Folder + 'ISI'+ '*J=100.' + SimuIdentifier1)
	PopRateI0_10 = glob.glob(Folder + 'PopRate'+ '*J=100.' + SimuIdentifier1)
	RasterI0_10 = glob.glob(Folder + 'Raster'+ '*J=100.' + SimuIdentifier1)
	CurrentI0_10 = glob.glob(Folder + 'FullMonitor'+ '*J=100.' + SimuIdentifier1)

	ISII0_19 = glob.glob(Folder + 'ISI'+ '*J=160.' + SimuIdentifier1)
	PopRateI0_19 = glob.glob(Folder + 'PopRate'+ '*J=160.' + SimuIdentifier1)
	RasterI0_19 = glob.glob(Folder + 'Raster'+ '*J=160.' + SimuIdentifier1)
	CurrentI0_19 = glob.glob(Folder + 'FullMonitor'+ '*J=160.' + SimuIdentifier1)
else:
	ISII0_10 = glob.glob(Folder + 'ISI'+ '*J=180.' + SimuIdentifier1)
	PopRateI0_10 = glob.glob(Folder + 'PopRate'+ '*J=180.' + SimuIdentifier1)
	RasterI0_10 = glob.glob(Folder + 'Raster'+ '*J=180.' + SimuIdentifier1)
	CurrentI0_10 = glob.glob(Folder + 'FullMonitor'+ '*J=180.' + SimuIdentifier1)

	ISII0_19 = glob.glob(Folder + 'ISI'+ '*J=180.' + SimuIdentifier1)
	PopRateI0_19 = glob.glob(Folder + 'PopRate'+ '*J=180.' + SimuIdentifier1)
	RasterI0_19 = glob.glob(Folder + 'Raster'+ '*J=180.' + SimuIdentifier1)
	CurrentI0_19 = glob.glob(Folder + 'FullMonitor'+ '*J=180.' + SimuIdentifier1)

ISIAux10 = np.loadtxt(ISII0_10[0],comments='#',skiprows=1)
t, PopAux10 = np.loadtxt(PopRateI0_10[0],comments='#',skiprows=1,unpack=True)
tSpikes10, IdxNeuAux10 = np.loadtxt(RasterI0_10[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux10, u_0Aux10 = np.loadtxt(CurrentI0_10[0],comments='#',skiprows=1,unpack=True)

ISIAux19 = np.loadtxt(ISII0_19[0],comments='#',skiprows=1)
t, PopAux19 = np.loadtxt(PopRateI0_19[0],comments='#',skiprows=1,unpack=True)
tSpikes19, IdxNeuAux19 = np.loadtxt(RasterI0_19[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux19, u_0Aux19 = np.loadtxt(CurrentI0_19[0],comments='#',skiprows=1,unpack=True)

Timestep=1
maxRateHyper=150
maxRateShunt=150
maxRate=maxRateHyper
xcoordinate=1905
tmin=1800
tmax=1900

NumbNeurons=3000.
Nu_0=30.

width = 7. # in inches
height = 10./2# in inches
Fig = plt.figure(figsize=[width,height]);

tmp = plt.subplot(2,4,1)
ax1 = plt.gca()
l = plt.hist(ISIAux10, bins=np.arange(0,200,10), histtype='stepfilled',range=(min(ISIAux10)*0.99,max(ISIAux10)*0.8), density=True, color='tab:blue')			
plt.title('A1',loc='left')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,2)
ax1 = plt.gca()
xlim_right=tSpikes10[-1]
xlim_left=tSpikes10[np.argmin(abs(tSpikes10-tmin))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes10, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes10)) )			
plt.title('A2',loc='left')
plt.ylabel('Counts')
plt.ylim([0,maxRate])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,3)
ax1 = plt.gca()
plt.scatter(tSpikes10, IdxNeuAux10, s=15.0, color ='black', marker = '|')
plt.title('A3',loc='left')
#plt.ylim([-0.5,501.5])
plt.ylim([-0.5,51.5])
plt.xlim([tmin,1900])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(2,4,4)
ax32 = plt.gca()
ax32.plot(t, v_0Aux10, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['right'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])
ax32.set_xlim([tmin,1900])
ax32.yaxis.label.set_color('blue')
ax32.tick_params(axis='y', colors='blue')

tmp = plt.subplot(2,4,5)
ax1 = plt.gca()
l = plt.hist(ISIAux19, bins=edgesbin, histtype='stepfilled',range=(min(ISIAux19)*0.99,max(ISIAux19)*1.01), density=True, color='tab:blue' )			
plt.title('B1',loc='left')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,6)
ax1 = plt.gca()
xlim_right=tSpikes19[-1]
xlim_left=tSpikes19[np.argmin(abs(tSpikes19-tmin))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes19, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes19)) )		
plt.title('B2',loc='left')
plt.ylabel('Counts')
plt.ylim([0,maxRate])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,7)
ax1 = plt.gca()
plt.scatter(tSpikes19, IdxNeuAux19, s=15.0, color ='black', marker = '|')
plt.title('B3',loc='left')
#plt.ylim([-0.5,501.5])
plt.ylim([-0.5,51.5])
plt.xlim([tmin,1900])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(2,4,8)
ax32 = plt.gca()
ax32.plot(t, v_0Aux19, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['right'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])
ax32.yaxis.label.set_color('blue')
ax32.tick_params(axis='y', colors='blue')
ax32.set_xlim([tmin,1900])

Fig.tight_layout()
Fig.subplots_adjust(wspace=0.7, hspace=0.3)
plt.savefig("FigureVia"+ Synapse + ".svg")
plt.savefig("FigureVia"+ Synapse + ".epsc")
