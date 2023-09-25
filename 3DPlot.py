import matplotlib.pyplot as plt
import numpy as np
import glob, os
from matplotlib.ticker import LinearLocator

sigmaAux = np.around(np.linspace(0.011,10.011,51, endpoint=True),4)
JAux = [1.,20.,40.,60.,80.,100.,120.,140.,160.,180.,250.,300.]
I0Aux= [30.,40.,45.,50.,55.,60.,65.,70.,75.,80.]

#FileIdentifier = "CriticalTransitionLine*IzhiTyII*0.0Hz.txt"
FileIdentifier = "CriticalTransitionLine*GuillemVarDriveHyper*.txt"

List=[]
for file1 in glob.glob(FileIdentifier):
	Matrix = np.loadtxt(file1);
	List.append(Matrix)
MatrixCritLine = np.row_stack(List)
np.savetxt("3DCriticalSurfGuillemHyper.txt",MatrixCritLine,fmt="%d %1.4f %d")
quit()

X = MatrixCritLine[:,1]
Y = MatrixCritLine[:,2]
Z = np.array( list(zip(*[iter(MatrixCritLine[:,0])]*1)) )

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# make the panes transparent
ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
# Customize the z axis.
ax.set_zlim(min(I0Aux)-0.1, max(I0Aux)+0.1)
ax.set_ylim(min(JAux)-0.1, max(JAux)+0.1)
#ax.zaxis.set_major_locator(LinearLocator(10))
# A StrMethodFormatter is used automatically
ax.zaxis.set_major_formatter('{x:.0f}')
ax.set_xlabel('Noise amplitude (mV)')
ax.set_ylabel('Synaptic strength (mV)')
ax.set_zlabel('External Input (mV)')

surf = ax.scatter(X, Y, Z,cmap='jet')
plt.savefig("IzhiTypeI.svg")

for ii in range(0,360,45):
	ax.view_init(elev=5., azim=ii)
	plt.savefig("IzhiTypeI%d.png" % ii)

plt.figure()
for i in range(0,len(glob.glob(FileIdentifier))):
	plt.plot(MatrixCritLine[i*6:5+i*6,1],MatrixCritLine[i*6:5+i*6,2],'x-',label='I0={I0}mV'.format(I0=I0Aux[i]))
plt.legend()
plt.savefig("IzhiTypeIPlaneJsigma.png")

plt.figure()
for i in range(0,len(JAux)):
	plt.plot(MatrixCritLine[i*8:8+i*8,1],MatrixCritLine[i*8:8+i*8,0],'x-',label='J={J}mV'.format(J=JAux[i]))
plt.ylim(min(I0Aux)-0.1, max(I0Aux)+0.1)
plt.legend()
plt.savefig("IzhiTypeIPlaneI0sigma.png")

