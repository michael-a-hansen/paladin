import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt


def load_matrix_market(matrixfile):
    """Load a matrix market file into a numpy array

    Keyword arguments:
    matrixfile -- the path to the matrix
    """
    return scio.mmread(source=matrixfile).toarray()


def load_spectrum(matrixfile):
    """Load a spectrum from a matrix file path and return a complex numpy array of eigenvalues.

    Keyword arguments:
    matrixfile -- the path to the matrix that was decomposed
    """
    data = np.loadtxt(fname=matrixfile + '.spectrum')
    return data[:, 0] + 1j * data[:, 1]


def load_leftvecs(matrixfile):
    """Load left eigenvectors from a matrix file path and return a complex numpy array of eigenvalues.

    Keyword arguments:
    matrixfile -- the path to the matrix that was decomposed
    """
    return load_matrix_market(matrixfile=matrixfile + '.leftvecs')


def load_rightvecs(matrixfile):
    """Load right eigenvectors from a matrix file path and return a complex numpy array of eigenvalues.

    Keyword arguments:
    matrixfile -- the path to the matrix that was decomposed
    """
    return load_matrix_market(matrixfile=matrixfile + '.rightvecs')


def plot_spectrum_on_complex_plane(matrixfile, axislimit, marker='ro'):
    """Load a spectrum from a matrix file path and plot its eigenvalues on the complex plane.

    Keyword arguments:
    matrixfile -- the path to the matrix that was decomposed
    axislimit -- the extent of the complex plane, for which axis limits will be [-axislimit, axislimit]^2
    marker -- the marker specification for the eigenvalues (default 'ro')

    Note that this loads the spectrum, draws and shows a matplotlib figure when called, so performance will be slow.
    If making a movie, this will NOT be optimal.
    Otherwise this should be fine.
    """
    w = load_spectrum(matrixfile=matrixfile)
    plt.plot(np.real(w), np.imag(w), marker)

    plt.plot([-axislimit, axislimit], [0, 0], 'k-')
    plt.plot([0, 0], [-axislimit, axislimit], 'k-')
    plt.xlim([-axislimit, axislimit])
    plt.ylim([-axislimit, axislimit])
    plt.grid(True)
    plt.xlabel('real')
    plt.ylabel('imag')
    plt.show()
