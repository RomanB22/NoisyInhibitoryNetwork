from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np

sigma = np.around(np.linspace(0.011,10.011,51, endpoint=True),4)
J = [1.,40.,80.,100.,120.,140.]
I0 = [0.,10.,20.,30.,40.,50.]

FileIdentifier = "3DCriticalSurfBisIzhiTyI.txt"
FileIdentifier2 = "3DCritSurfBis.txt"
FileIdentifier17Hz = "3DCriticalLine17Hz.txt"
FileIdentifier30Hz = "3DCriticalLine30Hz.txt"
FileIdentifier17HzTyI = "3DCriticalLine17HzTyI.txt"
FileIdentifier30HzTyI = "3DCriticalLine30HzTyI.txt"


CritLine17Hz = np.loadtxt(FileIdentifier17Hz)
CritLine30Hz = np.loadtxt(FileIdentifier30Hz)
CritLine17HzTyI = np.loadtxt(FileIdentifier17HzTyI)
CritLine30HzTyI = np.loadtxt(FileIdentifier30HzTyI)
MatrixCritLine = np.loadtxt(FileIdentifier)
MatrixCritLine2 = np.loadtxt(FileIdentifier2)

#CritLine17Hz = CritLine17Hz[CritLine17Hz[:,1]<12]
#CritLine30Hz = CritLine30Hz[CritLine30Hz[:,1]<12]
#MatrixCritLine2 = MatrixCritLine2[MatrixCritLine2[:,1]<12]

sigmaCrit = MatrixCritLine[:,1]
JAux = MatrixCritLine[:,2]
I0Aux = MatrixCritLine[:,0]

sigmaCrit2 = MatrixCritLine2[:,1]
JAux2 = MatrixCritLine2[:,2]
I0Aux2 = MatrixCritLine2[:,0]


# Plot the 3D surface
#I0Aux = np.ma.masked_where(I0Aux < 0.1, I0Aux)
#ax.plot_trisurf(sigmaCrit,JAux,I0Aux,edgecolor='royalblue', linewidth=0.2, alpha=0.3)
#X, Y = np.meshgrid(JAux,I0Aux)
#Z = np.zeros(np.shape(X))

#for i in range(np.shape(Z)[0]):
#	Z[i,i]=sigmaCrit[i]
#Z = np.ma.masked_where(Z < 0.1, Z)
#surf = ax.plot_surface(X, Y, Z, cmap='coolwarm', linewidth=0, antialiased=False)
#ax.plot_wireframe(X, Y, Z, rstride=2,  cstride=2,color='green')
#Z1 = np.array( list(zip(*[iter(I0Aux)]*1)) )
#surf = ax.scatter(sigmaCrit,JAux,Z1, cmap='jet')
#Z1 = np.array( list(zip(*[iter(sigmaCrit)]*1)) )
#surf = ax.scatter(JAux,I0Aux,Z1, cmap='jet')
#ax.set_zlim(0,30)

angleview=30
for ii in range(0,360,angleview):
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
	ax.view_init(elev=30., azim=ii)
	#ax.plot_trisurf(JAux,I0Aux,sigmaCrit,edgecolor='royalblue', linewidth=0.05, alpha=0.3)
	ax.plot_trisurf(JAux2,I0Aux2,sigmaCrit2,edgecolor='royalblue', linewidth=0.05, alpha=0.3)
	ax.plot(CritLine17Hz[:,2],CritLine17Hz[:,0],CritLine17Hz[:,1],"r",label=r'$\nu_0$=17 Hz Type II')
	ax.plot(CritLine30Hz[:,2],CritLine30Hz[:,0],CritLine30Hz[:,1],"g",label=r'$\nu_0$=30 Hz Type II')
	#ax.plot(CritLine17HzTyI[:,2],CritLine17HzTyI[:,0],CritLine17HzTyI[:,1],label=r'$\nu_0$=17 Hz Type I')
	#ax.plot(CritLine30HzTyI[:,2],CritLine30HzTyI[:,0],CritLine30HzTyI[:,1],label=r'$\nu_0$=30 Hz Type I')
	ax.text(75, 25, 3., 'Synchronous',(1,1,0))
	ax.text(5, 10, 15, 'Asynchronous',(1,1,0))
	ax.set_xlabel("Total synaptic strenght (mV)")
	ax.set_ylabel("Bias Current(mV)")
	ax.set_zlabel("Noise Intensity (mV)")
	plt.legend()
	plt.savefig("Transitions%d.svg" % ii)
plt.show()

