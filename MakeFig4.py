import numpy as np
import matplotlib.pyplot as plt
import glob, os
import matplotlib as mlb

def Reorder (tSpks, IdxSpks):
    LastSpk = []
    for i in np.unique(IdxSpks):
        IdxAux = np.argwhere(IdxSpks == i)
        LastSpk.append( [ i ,tSpks[ IdxAux[-1] ][0] ] )
    Mat = np.asarray( LastSpk )
    Idx = np.argsort(Mat[:,1]) 
    for j in Idx:
        OrderedtSpks = 1 
    return OrderedtSpks, OrderedIdxSpks 

mlb.rcParams['font.serif'] = "Arial"
mlb.rcParams['font.size'] = "10"

Folder = '../CSV_IzhiTyII/*'
SimuIdentifier1 = '*sigma=1.76*N=3000*nu0=30*'#0.2
SimuIdentifier2 = '*sigma=0.4*N=3000*nu0=30*'
SimuIdentifier3 = '*sigma=0.1*N=3000*nu0=30*'

ISII0_10 = glob.glob(Folder + 'ISI'+ '*J=1.6' + SimuIdentifier1)
PopRateI0_10 = glob.glob(Folder + 'PopRate'+ '*J=1.6' + SimuIdentifier1)
RasterI0_10 = glob.glob(Folder + 'Raster'+ '*J=1.6' + SimuIdentifier1)
CurrentI0_10 = glob.glob(Folder + 'FullMonitor'+ '*J=1.6' + SimuIdentifier1)

ISII0_19 = glob.glob(Folder + 'ISI'+ '*J=1.6' + SimuIdentifier2)
PopRateI0_19 = glob.glob(Folder + 'PopRate'+ '*J=1.6' + SimuIdentifier2)
RasterI0_19 = glob.glob(Folder + 'Raster'+ '*J=1.6' + SimuIdentifier2)
CurrentI0_19 = glob.glob(Folder + 'FullMonitor'+ '*J=1.6' + SimuIdentifier2)

ISII0_31 = glob.glob(Folder + 'ISI'+ '*J=1.6' + SimuIdentifier3)
PopRateI0_31 = glob.glob(Folder + 'PopRate'+ '*J=1.6' + SimuIdentifier3)
RasterI0_31 = glob.glob(Folder + 'Raster'+ '*J=1.6' + SimuIdentifier3)
CurrentI0_31 = glob.glob(Folder + 'FullMonitor'+ '*J=1.6' + SimuIdentifier3)

ISIAux10 = np.loadtxt(ISII0_10[0],comments='#',skiprows=1)
t, PopAux10 = np.loadtxt(PopRateI0_10[0],comments='#',skiprows=1,unpack=True)
tSpikes10, IdxNeuAux10 = np.loadtxt(RasterI0_10[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux10,v_1,v_2,v_3,v_4, u_0Aux10,u_1,u_2,u_3,u_4, Inotnoisy_0Aux10,InN_1,InN_2,InN_3,InN_4 = np.loadtxt(CurrentI0_10[0],comments='#',skiprows=1,unpack=True)

ISIAux19 = np.loadtxt(ISII0_19[0],comments='#',skiprows=1)
t, PopAux19 = np.loadtxt(PopRateI0_19[0],comments='#',skiprows=1,unpack=True)
tSpikes19, IdxNeuAux19 = np.loadtxt(RasterI0_19[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux19,v_1,v_2,v_3,v_4, u_0Aux19,u_1,u_2,u_3,u_4, Inotnoisy_0Aux19,InN_1,InN_2,InN_3,InN_4 = np.loadtxt(CurrentI0_19[0],comments='#',skiprows=1,unpack=True)

ISIAux31 = np.loadtxt(ISII0_31[0],comments='#',skiprows=1)
t, PopAux31 = np.loadtxt(PopRateI0_31[0],comments='#',skiprows=1,unpack=True)
tSpikes31, IdxNeuAux31 = np.loadtxt(RasterI0_31[0],comments='#',skiprows=1,unpack=True)	
t, v_0Aux31,v_1,v_2,v_3,v_4, u_0Aux31,u_1,u_2,u_3,u_4, Inotnoisy_0Aux31,InN_1,InN_2,InN_3,InN_4 = np.loadtxt(CurrentI0_31[0],comments='#',skiprows=1,unpack=True)

print(ISII0_10[0],ISII0_19[0],ISII0_31[0])

#print(tSpikes19, IdxNeuAux19)
#tSpikes19, IdxNeuAux19 = Reorder(tSpikes19,IdxNeuAux19)

Timestep=1
maxRate=400
xcoordinate=1805

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
height = 10.# in inches
Fig = plt.figure(figsize=[width,height]);

tmp = plt.subplot(3,4,1)
ax1 = plt.gca()
l = plt.hist(ISIAux10, bins=30, histtype='stepfilled',range=(min(ISIAux10)*0.99,max(ISIAux10)*1.01), density=True, color='tab:blue')			
plt.title('A1',loc='left')

plt.xlim([0,70])
ax1.spines['left'].set_visible(False)
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(3,4,2)
ax1 = plt.gca()
xlim_right=tSpikes10[-1]
xlim_left=tSpikes10[np.argmin(abs(tSpikes10-1700))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes10, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes10)) )			
plt.title('A2',loc='left')

plt.ylim([0,maxRate])
ax1.spines['left'].set_visible(False)
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(3,4,3)
ax1 = plt.gca()

plt.scatter(tSpikes10, IdxNeuAux10, s=15.0, color ='black', marker = '|')
plt.title('A3',loc='left')
plt.ylim([-0.5,70.5])
plt.xlim([1700,1800])
ax1.set_axis_off()

tmp = plt.subplot(3,4,4)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux10+3.59375, 'g-',label='Total')
ax3.plot(t, 3.59375 + np.mean(Inotnoisy_0Aux10) *np.ones(np.shape(t)), 'b-', linewidth=0.5,label='Mean')

ax3.plot(t, AH*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)

plt.title('A4',loc='left')
ax3.set_ylabel('Net Current (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.ylim(1,4.5)
ax3.spines['bottom'].set_visible(False)
ax3.set_xticks([])
ax3.set_xlim([1700,1800])
ax32 = ax3.twinx()
ax32.plot(t, v_0Aux10, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['left'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])



tmp = plt.subplot(3,4,5)
ax1 = plt.gca()
l = plt.hist(ISIAux19, bins=30, histtype='stepfilled',range=(min(ISIAux19)*0.99,max(ISIAux19)*1.01), density=True, color='tab:blue' )			
plt.title('B1',loc='left')

plt.xlim([0,70])
ax1.spines['left'].set_visible(False)
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(3,4,6)
ax1 = plt.gca()
xlim_right=tSpikes19[-1]
xlim_left=tSpikes19[np.argmin(abs(tSpikes19-1700))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes19, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue',weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes19)) )		
plt.title('B2',loc='left')

plt.ylim([0,maxRate])
ax1.spines['left'].set_visible(False)
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(3,4,7)
ax1 = plt.gca()
plt.scatter(tSpikes19, IdxNeuAux19, s=15.0, color ='black', marker = '|')
plt.title('B3',loc='left')
plt.ylim([-0.5,70.5])
plt.xlim([1700,1800])
ax1.set_axis_off()

tmp = plt.subplot(3,4,8)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux19+3.984375, 'g-',label='Total')
ax3.plot(t, 3.984375 + np.mean(Inotnoisy_0Aux19) *np.ones(np.shape(t)), 'b-', linewidth=0.5,label='Mean')

ax3.plot(t, AH*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)

plt.title('B4',loc='left')
ax3.set_ylabel('Net Current (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
plt.ylim(1,4.5)
ax3.spines['bottom'].set_visible(False)
ax3.set_xticks([])
ax3.set_xlim([1700,1800])
ax32 = ax3.twinx()
ax32.plot(t, v_0Aux19, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['left'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])



tmp = plt.subplot(3,4,9)
ax1 = plt.gca()
l = plt.hist(ISIAux31, bins=30, histtype='stepfilled',range=(min(ISIAux31)*0.99,max(ISIAux31)*1.01), density=True, color='tab:blue')			
plt.title('C1',loc='left')

plt.xlim([0,70])
ax1.spines['left'].set_visible(False)
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(3,4,10)
ax1 = plt.gca()
xlim_right=tSpikes31[-1]
xlim_left=tSpikes31[np.argmin(abs(tSpikes31-1700))]
Npoints = int((xlim_right-xlim_left)/Timestep) # Timestep 1 ms 
bin_edges = np.linspace(xlim_left, xlim_right,Npoints+1,endpoint=True)

l = plt.hist(tSpikes31, histtype='bar', bins=bin_edges, range=(xlim_left, xlim_right), density=False, color='tab:blue' ,weights=Nu_0**2/(NumbNeurons)*np.ones(np.shape(tSpikes31)) )			
plt.title('C2',loc='left')
plt.xlabel('t (ms)')

plt.ylim([0,maxRate])
ax1.spines['left'].set_visible(False)
ax1.set_yticks([])
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

tmp = plt.subplot(3,4,11)
ax1 = plt.gca()
plt.scatter(tSpikes31, IdxNeuAux31, s=15.0, color ='black', marker = '|')
plt.title('C3',loc='left')
plt.ylim([-0.5,70.5])
plt.xlim([1700,1800])
ax1.set_axis_off()

tmp = plt.subplot(3,4,12)
ax3 = plt.gca()
ax3.plot(t, Inotnoisy_0Aux31+3.5937500, 'g-',label='Total')
ax3.plot(t, 3.59375 + np.mean(Inotnoisy_0Aux31) *np.ones(np.shape(t)), 'b-', linewidth=0.5,label='Mean')

ax3.plot(t, AH*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateAH, "AH", fontsize=11)
ax3.plot(t, SNP*np.ones(np.shape(t)), 'r--')
ax3.text(xcoordinate, ycoordinateSNP, "SNP", fontsize=11)

plt.title('C4',loc='left')
ax3.set_ylabel('Net Current (mV)')
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.set_xticks([])
plt.ylim(1,4.5)
ax3.set_xlim([1700,1800])

ax32 = ax3.twinx()
ax32.plot(t, v_0Aux31, 'tab:blue')
ax32.set_ylabel('$V_m$ (mV)', color='blue')
ax32.spines['top'].set_visible(False)
ax32.spines['left'].set_visible(False)
ax32.spines['bottom'].set_visible(False)
ax32.set_ylim([-90,10])


Fig.tight_layout()
Fig.subplots_adjust(wspace=0.9, hspace=0.3)
plt.savefig("FigureCOtoCOA_30Hz.svg")

