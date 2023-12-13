"""receiver.py

Script with methods for the OTFS receiver.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np


def qam_sym_decision(y:np.ndarray, mod_size:int):
    """Method to decide QAM symbol.
    
    Parameters
    ----------
    y : np.ndarray
        Received symbols.
    mod_size : int
        Size of the QAM modulation.
    
    Returns
    -------
    rx_dig_sym : np.ndarray
        Received digitally modulated symbols.
    """

    if mod_size == 4:
        symbols = np.array((-1-1j, -1+1j, 1-1j, 1+1j))
    elif mod_size == 16:
        symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j, -1+1j,
                            -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j, 3-1j, 3+1j,
                            3+3j))
    n_rows, n_cols = y.shape
    rx_dig_sym = np.zeros((n_rows, n_cols), dtype=np.complex128)
    for i in range(n_rows):
        for j in range(n_cols):
            min_idx = np.argmin(np.abs(symbols - y[i, j]))
            rx_dig_sym[i, j] = symbols[min_idx]

    return rx_dig_sym