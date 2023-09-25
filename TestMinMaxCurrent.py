import numpy as np

Matrix = np.loadtxt("I0_AdExTyIIAH/MatrixI0JSigma_AdExTyIIAH17Hz_N800")

Current = Matrix[:,0]
print(min(Current),max(Current))
