#execfile("ImportingPackages.py")
from ImportingPackages import *
import pandas
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

nu0=17
NAux=[3000]
NTrucho=[3000]#NAux
ModelDescriptor = "AdExTyIIAH"


cols=np.arange(1,53)
skipRows=1


CAdEx = 150.

FileSynchronyIndex = "SPC{N}_{nu}_CSV_{Folder}.csv".format(N=NTrucho[0],nu=nu0,Folder=ModelDescriptor)
FileMaxInput = "MapMaxInputCurrent3000_"+str(nu0)+"_CSV_"+ModelDescriptor+".csv"
FileMeanInput = "MapMeanInputCurrent3000_"+str(nu0)+"_CSV_"+ModelDescriptor+".csv"
FileActivityMask = "ActivityMask_{nu}_CSV_{Folder}.csv".format(nu=nu0,Folder=ModelDescriptor)
MaxCurrent = np.loadtxt(FileMaxInput,delimiter=",",skiprows=skipRows,usecols=cols)
MeanCurrent = np.loadtxt(FileMeanInput,delimiter=",",skiprows=skipRows,usecols=cols)
ActivityMask = np.loadtxt(FileActivityMask,delimiter=",",skiprows=skipRows,usecols=cols)


J = np.loadtxt(FileSynchronyIndex,delimiter=",",max_rows=1)
sigma = np.loadtxt(FileSynchronyIndex,delimiter=",",skiprows=skipRows,usecols=0)
ChiInfinite = np.loadtxt(FileSynchronyIndex,delimiter=",",skiprows=skipRows,usecols=cols)

FileCriticalLine = "CriticalTransitionLineCSV_"+ModelDescriptor+"_"+str(nu0)+"Hz.txt"
CriticalLineAux = np.loadtxt(FileCriticalLine) 
sigmaCrit = CriticalLineAux[:,0]
JCrit = CriticalLineAux[:,1]
MaskSigma = np.zeros(np.shape(ChiInfinite))

for i in range(np.shape(JCrit)[0]):
	#print(JCrit[i],sigmaCrit[i],sigma[sigma < sigmaCrit[i]])	
	MaskSigma[:,i] = sigma < sigmaCrit[i]
MaskSigma[MaskSigma==0] = np.nan 

SynchronousPart = np.multiply(ChiInfinite,MaskSigma)

## Matplotlib Sample Code using 2D arrays via meshgrid
X, Y = np.meshgrid(sigma,J)
Z = SynchronousPart

ZMax = MaxCurrent
ZMean = MeanCurrent

#xnew, ynew = np.mgrid[min(J):max(J):0.2, min(sigma):max(sigma):(sigma[1]-sigma[0])/2]
#tck = interpolate.bisplrep(X, Y, Z, s=0.5)
#znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)

#znew = gaussian_filter(znew, sigma=1)

if ModelDescriptor == "IzhiTyI":
	levelsMaxCurr = [min(ZMax.flatten()),0.129321*10/0.09,max(ZMax.flatten())]
	levelsMeanCurr = [min(ZMean.flatten()),0.129321*10/0.09,max(ZMean.flatten())]
elif ModelDescriptor == "IzhiTyII":
	levelsMaxCurr = [min(ZMax.flatten()),0.1795*10.,0.2625*10.,0.4225*10.,max(ZMax.flatten())]
	levelsMeanCurr = [min(ZMean.flatten()),0.1795*10.,0.2625*10.,0.4225*10.,max(ZMean.flatten())]
elif ModelDescriptor == "AdExTyIIAH":
	factor=10./np.array(CAdEx)
	levelsMaxCurr = [min(ZMax.flatten()),factor*165,factor*184.8613,factor*185.2901,max(ZMax.flatten())]
	levelsMeanCurr = [min(ZMean.flatten()),factor*165,factor*184.8613,factor*185.2901,max(ZMean.flatten())]
elif ModelDescriptor == "LIF":
        levelsMaxCurr = [min(ZMax.flatten()),20.,max(ZMax.flatten())]
        levelsMeanCurr = [min(ZMean.flatten()),20.,max(ZMean.flatten())]





#xnew, ynew = np.mgrid[min(J):max(J):0.2, min(sigma):max(sigma):(sigma[1]-sigma[0])/2]
#tck = interpolate.bisplrep(X, Y, Z, s=0.5)
#znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)

#znew = gaussian_filter(znew, sigma=1)

fig = plt.figure()

ax = fig.gca()
#cont = ax.contourf(np.log(X), np.log(Y), np.transpose(Z), levels=np.linspace(0.15,1.,40,endpoint=True), cmap=cm.viridis)

cont = ax.contourf(np.log(X), np.log(Y), np.transpose(Z), levels=np.linspace(0.15,1.,40,endpoint=True), cmap=cm.viridis)


cbar = fig.colorbar(cont, shrink=0.7, aspect=10)

cont1 = ax.contour(np.log(X), np.log(Y), np.transpose(ZMax), levels=levelsMaxCurr, linestyles='dashdot')
cont2 = ax.contour(np.log(X), np.log(Y), np.transpose(ZMean), levels=levelsMeanCurr)
niveles=[0.1,2970,3001]#np.linspace(min(ActivityMask[ActivityMask>0.]),NAux[0],3,endpoint=True)

cont3 = ax.contour(np.log(X), np.log(Y), np.transpose(ActivityMask), levels=[0.1,2970,3001],linestyles='dashdot')

ax.plot(np.log(sigmaCrit),np.log(JCrit),"rx-")

ax.set_xticks([np.log(0.01),np.log(0.02),np.log(0.03),np.log(0.04),np.log(0.05),np.log(0.06),np.log(0.07),np.log(0.08),np.log(0.09),np.log(0.1),np.log(0.2),np.log(0.3),np.log(0.4),np.log(0.5),np.log(0.6),np.log(0.7),np.log(0.8),np.log(0.9),np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10)])


ax.set_yticks([np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10),np.log(20),np.log(30),np.log(40),np.log(50),np.log(60),np.log(70),np.log(80),np.log(90),np.log(100),np.log(110),np.log(120),np.log(130),np.log(140),np.log(153)])

ax.set_xticklabels(["0.01","0.02","","","","","","","","0.1","0.2","","","","","","","","1","2","","","","","","","","10"])
ax.set_yticklabels(["1","2","3","4","5","6","","","","10","20","30","40","50","60","","","","100","","","","","153"])

ax.set_xlim(np.log(0.01), np.log(max(sigma)))

ax.set_ylim(np.log(min(J)), np.log(max(J)))

ax.set_xlabel("Noise amplitude (mV)")
ax.set_ylabel("Total synatic strength (mV)")
#cbar.set_clim(0, 1)
plt.tight_layout()

plt.savefig("TotalMap{N}_{nu}_{Folder}.svg".format(N=3000,nu=nu0,Folder=ModelDescriptor))
plt.savefig("TotalMap{N}_{nu}_{Folder}.eps".format(N=3000,nu=nu0,Folder=ModelDescriptor))
