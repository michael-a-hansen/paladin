import pypp.pypp as pp

# path to the matrix that was decomposed
matrixfile = 'example-data/matrix'

# print the eigenvalues
w = pp.load_spectrum(matrixfile=matrixfile)
print('eigenvalues: ')
print(w)

# plot the eigenvalues on the complex plane
pp.plot_spectrum_on_complex_plane(matrixfile=matrixfile, axislimit=2.1, marker='ro')
