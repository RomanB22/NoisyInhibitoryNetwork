#execfile("ImportingPackages.py")
from ImportingPackages import *
import pandas
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp1d

nu0=30
NAux=[3000]
ModelDescriptor = "IzhiTyI"
CAdEx = 150.


cols=np.arange(1,53)

FileSynchronyIndex = "SPC{N}_{nu}_CSV_{Folder}.csv".format(N=NAux[0],nu=nu0,Folder=ModelDescriptor)
J = np.loadtxt(FileSynchronyIndex,delimiter=",",max_rows=1)
sigma = np.loadtxt(FileSynchronyIndex,delimiter=",",skiprows=1,usecols=0)

FileExtInput = "ExtCurrent3000_"+str(nu0)+"_CSV_"+ModelDescriptor+".csv"
ExtCurrent = np.loadtxt(FileExtInput,delimiter=",",skiprows=1,usecols=cols)

FileCriticalLine = "CriticalTransitionLineCSV_"+ModelDescriptor+"_"+str(nu0)+"Hz.txt"
CriticalLineAux = np.loadtxt(FileCriticalLine) 
sigmaCrit = CriticalLineAux[:,0]
JCrit = CriticalLineAux[:,1]
MaskSigma = np.zeros(np.shape(ExtCurrent))
MaskSigma = ExtCurrent >= -10.
#MaskSigma[MaskSigma==0] = 0 

SynchronousPart = np.multiply(ExtCurrent,MaskSigma)
## Matplotlib Sample Code using 2D arrays via meshgrid
X, Y = np.meshgrid(sigma,J)
#Z = ExtCurrent[ExtCurrent>=0.]
Z=SynchronousPart

levels=40

fig = plt.figure()
ax = fig.gca()
cont = ax.contourf(np.log(X), np.log(Y), np.transpose(Z), levels=levels, cmap=cm.viridis)

cbar = fig.colorbar(cont, shrink=0.7, aspect=10)
ax.set_xticks([np.log(0.01),np.log(0.02),np.log(0.03),np.log(0.04),np.log(0.05),np.log(0.06),np.log(0.07),np.log(0.08),np.log(0.09),np.log(0.1),np.log(0.2),np.log(0.3),np.log(0.4),np.log(0.5),np.log(0.6),np.log(0.7),np.log(0.8),np.log(0.9),np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10)])
ax.set_yticks([np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10),np.log(20),np.log(30),np.log(40),np.log(50),np.log(60),np.log(70),np.log(80),np.log(90),np.log(100),np.log(110),np.log(120),np.log(130),np.log(140),np.log(153)])

ax.plot(np.log(sigmaCrit),np.log(JCrit),"rx-")

ax.set_xticklabels(["0.01","0.02","","","","","","","","0.1","0.2","","","","","","","","1","2","","","","","","","","10"])
ax.set_yticklabels(["1","2","3","4","5","6","","","","10","20","30","40","50","60","","","","100","","","","","153"])

ax.set_xlim(np.log(0.1), np.log(max(sigma)))
ax.set_ylim(np.log(min(J)), np.log(max(J)))


#JAux = [1.3,6.3,20.,39.8,100.]
#sigmaAux = [0.1122,0.261,0.5309,1.011,1.511,2.511,4.011,4.217,5.511,8.761,10.011]

sigmaAux = np.sort(np.concatenate((np.around(np.logspace(-2,1,41,endpoint=True),4),np.linspace(0.011,10.011,41, endpoint=True))))
JAux = np.sort(np.concatenate((np.around(np.logspace(0,2,21,endpoint=True),1),np.linspace(3,153,31, endpoint=True))))

#JAuxLines = np.log(JAux)
#sigmaAuxLines = np.log(sigmaAux)

#X2, Y2 = np.meshgrid(sigmaAuxLines,JAuxLines)
#plt.plot(X2,Y2,'cd')

ax.set_xlabel("Noise amplitude (mV)")
ax.set_ylabel("Total synatic strength (mV)")
#cbar.set_clim(min(Z.flatten()),max(Z.flatten()))
plt.tight_layout()

plt.savefig("ExtInput{N}_{nu}_{Folder}.svg".format(N=3000,nu=nu0,Folder=ModelDescriptor))

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})



# I want to interpolate over Sigma for each J
Icrit =[]
for i in range(np.size(JAux)):
	f2 = interp1d(sigmaAux, SynchronousPart[:,i], kind='cubic')
	if sigmaCrit[i] <= max(sigmaAux):
		Icrit.append(f2(sigmaCrit[i]))
	else:	
		Icrit.append(0)

Icrit = np.array(Icrit)
Icrit[Icrit<-10]=0

ax.plot(np.log(sigmaCrit),np.log(JCrit),Icrit,"rx-")

surf = ax.plot_surface(np.log(X), np.log(Y), np.transpose(Z), cmap=cm.viridis, linewidth=0, antialiased=False)

cbar = fig.colorbar(cont, shrink=0.7, aspect=10)
ax.set_xticks([np.log(0.01),np.log(0.02),np.log(0.03),np.log(0.04),np.log(0.05),np.log(0.06),np.log(0.07),np.log(0.08),np.log(0.09),np.log(0.1),np.log(0.2),np.log(0.3),np.log(0.4),np.log(0.5),np.log(0.6),np.log(0.7),np.log(0.8),np.log(0.9),np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10)])
ax.set_yticks([np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10),np.log(20),np.log(30),np.log(40),np.log(50),np.log(60),np.log(70),np.log(80),np.log(90),np.log(100),np.log(110),np.log(120),np.log(130),np.log(140),np.log(153)])

ax.set_xticklabels(["0.01","","","","","","","","","0.1","","","","","","","","","1","","","","","","","","","10"])
ax.set_yticklabels(["1","","","","","","","","","10","","","","","","","","","100","","","","","153"])

ax.set_xlim(np.log(0.1), np.log(max(sigma)))
ax.set_ylim(np.log(min(J)), np.log(max(J)))

ax.set_xlabel("Noise")
ax.set_ylabel("Coupling")

ax.view_init(azim=-100, elev=30)

plt.savefig("3D_ExtInput{N}_{nu}_{Folder}.svg".format(N=3000,nu=nu0,Folder=ModelDescriptor))
