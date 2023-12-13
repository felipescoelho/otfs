"""utils.py

Script with utilitary functions.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np


def dft_mat(N):
    """Method to generate the DFT matrix.
    
    Parameters
    ----------
    N : int
        Number of bins in the DFT matrix.
    
    Returns
    -------
    dft_mtx : np.ndarray
        The DFT matrix.
    """

    dft_mtx = np.zeros((N, N), dtype=np.complex128)
    for k in range(N):
        for l in range(N):
            dft_mtx[k, l] = np.exp(-1j*2*np.pi*k*l/N)/np.sqrt(N)

    return dft_mtx