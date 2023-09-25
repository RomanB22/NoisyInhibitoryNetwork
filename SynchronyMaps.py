#execfile("ImportingPackages.py")
from ImportingPackages import *
import pandas
from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

nu0=30
ModelDescriptor = "IzhiTyI"


cols=np.arange(1,53)

FileSynchronyIndex = "HeatMapChiIfty_"+str(nu0)+"_CSV_"+ModelDescriptor+".csv"

J = np.loadtxt(FileSynchronyIndex,delimiter=",",max_rows=1)
sigma = np.loadtxt(FileSynchronyIndex,delimiter=",",skiprows=1,usecols=0)
ChiInfinite = np.loadtxt(FileSynchronyIndex,delimiter=",",skiprows=1,usecols=cols)



FileCriticalLine = "CriticalTransitionLineCSV_"+ModelDescriptor+"_"+str(nu0)+"Hz.txt"
CriticalLineAux = np.loadtxt(FileCriticalLine) 
sigmaCrit = CriticalLineAux[:,0]
JCrit = CriticalLineAux[:,1]
MaskSigma = np.zeros(np.shape(ChiInfinite))
for i in range(np.shape(JCrit)[0]):	
		MaskSigma[:,i] = sigma < sigmaCrit[i]
MaskSigma[MaskSigma==0] = np.nan 

SynchronousPart = np.multiply(ChiInfinite,MaskSigma)

## Matplotlib Sample Code using 2D arrays via meshgrid
X, Y = np.meshgrid(sigma,J)
Z = SynchronousPart

#xnew, ynew = np.mgrid[min(J):max(J):0.2, min(sigma):max(sigma):(sigma[1]-sigma[0])/2]
#tck = interpolate.bisplrep(X, Y, Z, s=0.5)
#znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)

#znew = gaussian_filter(znew, sigma=1)

fig = plt.figure()

ax = fig.gca()
cont = ax.contourf(np.log(X), np.log(Y), np.transpose(Z), levels=40, cmap=cm.viridis)

cbar = fig.colorbar(cont, shrink=0.7, aspect=10)
ax.plot(np.log(sigmaCrit),np.log(JCrit),"rx-")

ax.set_xticks([np.log(0.01),np.log(0.02),np.log(0.03),np.log(0.04),np.log(0.05),np.log(0.06),np.log(0.07),np.log(0.08),np.log(0.09),np.log(0.1),np.log(0.2),np.log(0.3),np.log(0.4),np.log(0.5),np.log(0.6),np.log(0.7),np.log(0.8),np.log(0.9),np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10)])


ax.set_yticks([np.log(1),np.log(2),np.log(3),np.log(4),np.log(5),np.log(6),np.log(7),np.log(8),np.log(9),np.log(10),np.log(20),np.log(30),np.log(40),np.log(50),np.log(60),np.log(70),np.log(80),np.log(90),np.log(100),np.log(110),np.log(120),np.log(130),np.log(140),np.log(153)])

ax.set_xticklabels(["0.01","0.02","","","","","","","","0.1","0.2","","","","","","","","1","2","","","","","","","","10"])
ax.set_yticklabels(["1","2","3","4","5","6","","","","10","20","30","40","50","60","","","","100","","","","","153"])

ax.set_xlim(np.log(0.1), np.log(max(sigma)))
ax.set_ylim(np.log(min(J)), np.log(max(J)))

ax.set_xlabel("Noise amplitude (mV)")
ax.set_ylabel("Total synatic strength (mV)")

cbar.set_clim(0, 1)
plt.tight_layout()

plt.savefig("PhaseDiagram{N}_{nu}_{Folder}.svg".format(N=3000,nu=nu0,Folder=ModelDescriptor))
