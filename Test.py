import numpy as np


params = np.loadtxt("params.dat"); # Parameters: EL, Rin and tauv, gNaf, rat (=gKf/gNaf), thm1, thh2, thn1, and gKv1f, tha1
lento = np.size(params,1)
print(lento)

A=np.zeros((lento,1))

for ChosenNeuronModelGuillem in [56,57]:#range(0,lento):
	ELs, Rins, tauvs, gNafs, rats, thm1s, thh2s, thn1s, gKv1fs, tha1s = params[:,ChosenNeuronModelGuillem]
	A[ChosenNeuronModelGuillem,:]=tauvs
	print(tauvs)
print("Max",max(A),np.argmax(A))
