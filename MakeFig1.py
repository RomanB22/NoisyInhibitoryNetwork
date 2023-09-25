import numpy as np
import matplotlib.pyplot as plt
import glob, os
import matplotlib as mlb

mlb.rcParams['font.serif'] = "Arial"
mlb.rcParams['font.size'] = "6"

Folder = '../CSV_IzhiTyIIVarDrive/*'
SimuIdentifier1 = '*J=79.4*sigma=10.0110*'
SimuIdentifier2 = '*N=3000*'

ISII0_10 = glob.glob(Folder + 'ISI'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=4.*-50.0*')
PopRateI0_10 = glob.glob(Folder + 'PopRate'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=4.*-50.0*')
RasterI0_10 = glob.glob(Folder + 'Raster'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=4.*-50.0*')
CurrentI0_10 = glob.glob(Folder + 'FullMonitor'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=4.*-50.0*')

ISII0_19 = glob.glob(Folder + 'ISI'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=10.*-50.0*')
PopRateI0_19 = glob.glob(Folder + 'PopRate'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=10.*-50.0*')
RasterI0_19 = glob.glob(Folder + 'Raster'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=10.*-50.0*')
CurrentI0_19 = glob.glob(Folder + 'FullMonitor'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=10.*-50.0*')

ISII0_31 = glob.glob(Folder + 'ISI'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=16.*-50.0*')
PopRateI0_31 = glob.glob(Folder + 'PopRate'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=16.*-50.0*')
RasterI0_31 = glob.glob(Folder + 'Raster'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=16.*-50.0*')
CurrentI0_31 = glob.glob(Folder + 'FullMonitor'+ SimuIdentifier1 + SimuIdentifier2 + 'I0=16.*-50.0*')

ISIAux10 = np.loadtxt(ISII0_10[0],comments='#',skiprows=1)
t, PopAux10 = np.loadtxt(PopRateI0_10[0],comments='#',skiprows=1,unpack=True)
tSpikes10, IdxNeuAux10 = np.loadtxt(RasterI0_10[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux10, u_0Aux10, Inotnoisy_0Aux10 = np.loadtxt(CurrentI0_10[0],comments='#',skiprows=1,unpack=True)

print(min(t),max(t),len(t))

ISIAux19 = np.loadtxt(ISII0_19[0],comments='#',skiprows=1)
t, PopAux19 = np.loadtxt(PopRateI0_19[0],comments='#',skiprows=1,unpack=True)
tSpikes19, IdxNeuAux19 = np.loadtxt(RasterI0_19[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux19, u_0Aux19, Inotnoisy_0Aux19 = np.loadtxt(CurrentI0_19[0],comments='#',skiprows=1,unpack=True)

print(min(t),max(t),len(t))

ISIAux31 = np.loadtxt(ISII0_31[0],comments='#',skiprows=1)
t, PopAux31 = np.loadtxt(PopRateI0_31[0],comments='#',skiprows=1,unpack=True)
tSpikes31, IdxNeuAux31 = np.loadtxt(RasterI0_31[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux31, u_0Aux31, Inotnoisy_0Aux31 = np.loadtxt(CurrentI0_31[0],comments='#',skiprows=1,unpack=True)

print(min(t),max(t),len(t));

Timestep=1
maxRate=100
xcoordinate=1605


factor=10 # Izhi Type II
SN = factor*0.4225
ycoordinateSN=SN
AH = factor*0.2625
ycoordinateAH=AH
SNP = factor*0.1795
ycoordinateSNP=SNP

NumbNeurons=3000.
Nu_0=30.

#width = 7. # in inches
#height = 10.# in inches
#Fig = plt.figure(figsize=[width,height])
Fig = plt.figure()

tmp = plt.subplot(3,4,1)
ax1 = plt.gca()
l = plt.hist(ISIAux10, bins=20, histtype='stepfilled',range=(0,200), density=True, color='tab:blue')			
plt.title('A1',loc='left')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,2)
ax1 = plt.gca()
xlim_right=tSpikes10[-1]
xlim_left=tSpikes10[np.argmin(abs(tSpikes10-1500))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes10, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes10)) )			
plt.title('A2',loc='left')
#plt.ylabel('Counts')
plt.ylim([0,maxRate])
ax1.spines['top'].set_visible(False)
ax1.spines['left'].set_visible(False)
ax1.spines['right'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,3)
ax1 = plt.gca()
plt.scatter(tSpikes10, IdxNeuAux10, s=5.0, color ='k', marker = '|')
plt.title('A3',loc='left')
#plt.ylim([-0.5,21.5])
plt.ylim([-0.5,70.5])
plt.xlim([1500,1600])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,4)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux10+4, 'g',label='Total')
ax3.plot(t, 4 + np.mean(Inotnoisy_0Aux10) *np.ones(np.shape(t)), 'b', linewidth=0.5,label='Mean')
#ax3.plot(t, SN*np.ones(np.shape(t)), 'r--')
#ax3.text(xcoordinate, ycoordinateSN, "SN", fontsize=11)
ax3.plot(t, AH*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)
#plt.legend()
plt.title('A4',loc='left')
ax3.set_ylabel('Total Synaptic Input (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.ylim([-25,5])
ax3.set_xlim([1500,1600])
ax32 = ax3.twinx()
ax32.plot(t, v_0Aux10, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['left'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])
ax32.spines['right'].set_color('blue')
ax32.yaxis.label.set_color('blue')
ax32.tick_params(axis='y', colors='blue')

tmp = plt.subplot(3,4,5)
ax1 = plt.gca()
l = plt.hist(ISIAux19, bins=20, histtype='stepfilled',range=(0,200), density=True, color='tab:blue' )			
plt.title('B1',loc='left')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,6)
ax1 = plt.gca()
xlim_right=tSpikes19[-1]
xlim_left=tSpikes19[np.argmin(abs(tSpikes19-1500))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes19, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes19)) )		
plt.title('B2',loc='left')
plt.ylabel('Counts')
plt.ylim([0,maxRate])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,7)
ax1 = plt.gca()
plt.scatter(tSpikes19, IdxNeuAux19, s=5.0, color ='k', marker = '|')
plt.title('B3',loc='left')
#plt.ylim([-0.5,21.5])
plt.ylim([-0.5,70.5])
plt.xlim([1500,1600])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,8)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux19+10, 'g',label='Total')
ax3.plot(t, 10 + np.mean(Inotnoisy_0Aux19) *np.ones(np.shape(t)), 'b', linewidth=0.5,label='Mean')
#ax3.plot(t, SN*np.ones(np.shape(t)), 'r--')
#ax3.text(xcoordinate, ycoordinateSN, "SN", fontsize=11)
ax3.plot(t, AH*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)
#plt.legend()
plt.title('B4',loc='left')
ax3.set_ylabel('Total Synaptic Input (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.ylim([-25,5])
ax3.set_xlim([1500,1600])
ax32 = ax3.twinx()
ax32.plot(t, v_0Aux19, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['left'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])
ax32.spines['right'].set_color('blue')
ax32.yaxis.label.set_color('blue')
ax32.tick_params(axis='y', colors='blue')

tmp = plt.subplot(3,4,9)
ax1 = plt.gca()
l = plt.hist(ISIAux31, bins=20, histtype='stepfilled',range=(0,200), density=True, color='tab:blue')			
plt.title('C1',loc='left')
plt.xlabel('ISI (ms)')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,10)
ax1 = plt.gca()
xlim_right=tSpikes31[-1]
xlim_left=tSpikes31[np.argmin(abs(tSpikes31-1500))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)
l = plt.hist(tSpikes31, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue' ,weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes31)) )			
plt.title('C2',loc='left')
plt.xlabel('t (ms)')
plt.ylabel('Counts')
plt.ylim([0,maxRate])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,11)
ax1 = plt.gca()
plt.scatter(tSpikes31, IdxNeuAux31, s=5.0, color ='black', marker = '|')
plt.title('C3',loc='left')
plt.xlabel('t (ms)')
#plt.ylim([-0.5,21.5])
plt.ylim([-0.5,70.5])
plt.xlim([1500,1600])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(3,4,12)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux31+13, 'g',label='Total')
ax3.plot(t, 13 + np.mean(Inotnoisy_0Aux31) *np.ones(np.shape(t)), 'b', linewidth=0.5,label='Mean')
#ax3.plot(t, SN*np.ones(np.shape(t)), 'r--')
#ax3.text(xcoordinate, ycoordinateSN, "SN", fontsize=11)
ax3.plot(t, AH*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)
#plt.legend()
plt.title('C4',loc='left')
ax3.set_ylabel('Total Synaptic Input (mV)')
ax3.set_xlabel('t (ms)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.ylim([-25,5])
ax3.set_xlim([1500,1600])
ax32 = ax3.twinx()
ax32.plot(t, v_0Aux31, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['left'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])
ax32.spines['right'].set_color('blue')
ax32.yaxis.label.set_color('blue')
ax32.tick_params(axis='y', colors='blue')


Fig.tight_layout()
#Fig.subplots_adjust(wspace=0.7, hspace=0.2)
plt.savefig("FigureSPAtoSPO.svg")
plt.savefig("FigureSPAtoSPO.eps")
