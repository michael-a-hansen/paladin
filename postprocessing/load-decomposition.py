import pypp.pypp as pp
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

# path to the matrix that was decomposed
matrixfile = 'example-data/matrix'

# load the matrix
matrix = pp.load_matrix_market(matrixfile=matrixfile)

# get the decomposition from paladin
w = pp.load_spectrum(matrixfile=matrixfile)
vl = pp.load_leftvecs(matrixfile=matrixfile)
vr = pp.load_rightvecs(matrixfile=matrixfile)

# get the decomposition from scipy.linalg
w2, vl2, vr2 = la.eig(a=matrix, left=True, right=True)

# print the dominant eigenpair
domidx = np.argmax(np.abs(w))
domidx2 = np.argmax(np.abs(w2))
print('paladin: ', w[domidx], ': ', vl[:, domidx])
print('scipy.linalg: ', w2[domidx2], ': ', vl2[:, domidx2])

# plot the eigenvectors in the complex plane
for i in np.arange(0, 20):
    plt.plot(np.real(vl[:, i]), np.imag(vl[:, i]), 'b')
    plt.plot(np.real(vl2[:, i]), np.imag(vl2[:, i]), 'r--')

plt.grid(True)
plt.xlabel('real')
plt.ylabel('imag')
plt.show()
