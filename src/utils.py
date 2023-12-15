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


def row_col_permutation_mat(num_doppler_bins:int, num_delay_bins:int):
    """Method to generate the row-column permutation matrix.
    
    Parameters
    ----------
    num_doppler_bins : int
        Number of Doppler bins.
    num_delay_bins : int
        Number of delay bins.
    
    Returns
    -------
    P_mat : np.ndarray
        Row-column permutation matrix.
    """

    num_sym_per_frame = int(num_doppler_bins*num_delay_bins)
    P_mat = np.zeros((num_sym_per_frame, num_sym_per_frame))
    for j in range(num_doppler_bins):
        for i in range(num_delay_bins):
            aux_mat = np.zeros((num_delay_bins, num_doppler_bins))
            aux_mat[i, j] = 1
            P_mat[int(j*num_delay_bins):int((j+1)*num_delay_bins),
                  int(i*num_doppler_bins):int((i+1)*num_doppler_bins)] = aux_mat
    
    return P_mat