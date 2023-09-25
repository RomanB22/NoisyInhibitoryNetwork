import numpy as np
import matplotlib.pyplot as plt
import glob, os
import matplotlib as mlb

mlb.rcParams['font.serif'] = "Arial"
mlb.rcParams['font.size'] = "10"

Folder = '../CSV_IzhiTyII/*'
SimuIdentifier1 = '*sigma=1.0110*N=3000*nu0=17*'
SimuIdentifier2 = '*sigma=0.1122*N=3000*nu0=17*'

ISII0_10 = glob.glob(Folder + 'ISI'+ '*J=2.0' + SimuIdentifier1)
PopRateI0_10 = glob.glob(Folder + 'PopRate'+ '*J=2.0' + SimuIdentifier1)
RasterI0_10 = glob.glob(Folder + 'Raster'+ '*J=2.0' + SimuIdentifier1)
CurrentI0_10 = glob.glob(Folder + 'FullMonitor'+ '*J=2.0' + SimuIdentifier1)

ISII0_19 = glob.glob(Folder + 'ISI'+ '*J=2.0' + SimuIdentifier2)
PopRateI0_19 = glob.glob(Folder + 'PopRate'+ '*J=2.0' + SimuIdentifier2)
RasterI0_19 = glob.glob(Folder + 'Raster'+ '*J=2.0' + SimuIdentifier2)
CurrentI0_19 = glob.glob(Folder + 'FullMonitor'+ '*J=2.0' + SimuIdentifier2)

ISIAux10 = np.loadtxt(ISII0_10[0],comments='#',skiprows=1)
t, PopAux10 = np.loadtxt(PopRateI0_10[0],comments='#',skiprows=1,unpack=True)
tSpikes10, IdxNeuAux10 = np.loadtxt(RasterI0_10[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux10,v_1,v_2,v_3,v_4, u_0Aux10,u_1,u_2,u_3,u_4, Inotnoisy_0Aux10,InN_1,InN_2,InN_3,InN_4 = np.loadtxt(CurrentI0_10[0],comments='#',skiprows=1,unpack=True)

ISIAux19 = np.loadtxt(ISII0_19[0],comments='#',skiprows=1)
t, PopAux19 = np.loadtxt(PopRateI0_19[0],comments='#',skiprows=1,unpack=True)
tSpikes19, IdxNeuAux19 = np.loadtxt(RasterI0_19[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux19,v_1,v_2,v_3,v_4, u_0Aux19,u_1,u_2,u_3,u_4, Inotnoisy_0Aux19,InN_1,InN_2,InN_3,InN_4 = np.loadtxt(CurrentI0_19[0],comments='#',skiprows=1,unpack=True)

Timestep=1
maxRate=700
xcoordinate=1805
tmin=1500

factor=10 # Izhi Type II
SN = factor*0.4225
ycoordinateSN=SN
AH = factor*0.2625
ycoordinateAH=AH
SNP = factor*0.1795
ycoordinateSNP=SNP

NumbNeurons=3000.
Nu_0=30.

width = 7. # in inches
height = 10./2# in inches
Fig = plt.figure(figsize=[width,height]);

tmp = plt.subplot(2,4,1)
ax1 = plt.gca()
l = plt.hist(ISIAux10, bins=np.arange(20,80,1), histtype='stepfilled',range=(min(ISIAux10)*0.99,max(ISIAux10)*0.8), density=True, color='tab:blue')			
plt.title('A1-J=7.9 mV,$\sigma$=1 mV',loc='left')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,2)
ax1 = plt.gca()
xlim_right=tSpikes10[-1]
xlim_left=tSpikes10[np.argmin(abs(tSpikes10-tmin))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes10, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='k',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes10)) )			
plt.title('A2',loc='left')
plt.ylabel('Counts')
plt.ylim([0,maxRate/5])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,3)
ax1 = plt.gca()
plt.scatter(tSpikes10, IdxNeuAux10, s=15.0, color ='black', marker = '|')
plt.title('A3',loc='left')
plt.ylim([21.5,42.5])
plt.xlim([tmin,1800])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(2,4,4)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux10+4.707, 'k-',label='Total')
ax3.plot(t, 4.707 + np.mean(Inotnoisy_0Aux10) *np.ones(np.shape(t)), 'k--', linewidth=0.5,label='Mean')
ax3.plot(t, SN*np.ones(np.shape(t)), 'k-.')
ax3.text(xcoordinate, ycoordinateSN, "SN", fontsize=11)
ax3.plot(t, AH*np.ones(np.shape(t)), 'k-.')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'k-.')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)
plt.legend()
plt.title('A4',loc='left')
ax3.set_ylabel('Net Current (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.ylim(-1,8)
ax3.set_xlim([tmin,1800])


tmp = plt.subplot(2,4,5)
ax1 = plt.gca()
l = plt.hist(ISIAux19, bins=30, histtype='stepfilled',range=(min(ISIAux19)*0.99,max(ISIAux19)*1.01), density=True, color='k' )			
plt.title('B1-J=148 mV,$\sigma$=0.1 mV',loc='left')
plt.ylabel('Counts')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,6)
ax1 = plt.gca()
xlim_right=tSpikes19[-1]
xlim_left=tSpikes19[np.argmin(abs(tSpikes19-tmin))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes19, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='k',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes19)) )		
plt.title('B2',loc='left')
plt.ylabel('Counts')
plt.ylim([0,maxRate])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(2,4,7)
ax1 = plt.gca()
plt.scatter(tSpikes19, IdxNeuAux19, s=15.0, color ='black', marker = '|')
plt.title('B3',loc='left')
plt.ylim([-0.5,21.5])
plt.xlim([tmin,1800])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)
plt.yticks([])

tmp = plt.subplot(2,4,8)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux19+2.8125, 'k-',label='Total')
ax3.plot(t, 2.8125 + np.mean(Inotnoisy_0Aux19) *np.ones(np.shape(t)), 'k--', linewidth=0.5,label='Mean')
ax3.plot(t, SN*np.ones(np.shape(t)), 'k-.')
ax3.text(xcoordinate, ycoordinateSN, "SN", fontsize=11)
ax3.plot(t, AH*np.ones(np.shape(t)), 'k-.')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'k-.')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)
plt.legend()
plt.title('B4',loc='left')
ax3.set_ylabel('Net Current (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.ylim(-220,5)
ax3.set_xlim([tmin,1800])

Fig.tight_layout()
Fig.subplots_adjust(wspace=0.7, hspace=0.3)
plt.savefig("Figure617Hz.svg")
