"""transmitter.py

Script with methods for the OTFS transmitter.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np


def qam_sym(num_sym:int, mod_size:int):
    """Method to generate random QAM symbols.
    
    Parameters
    ----------
    num_sym : int
        Number of symbols.
    mod_size : int
        Size of the QAM constellation.
    
    Returns
    -------
    signal_digmod : np.ndarray
        Digitally modulated signal.
    """

    if mod_size == 4:
        symbols = np.array((-1-1j, -1+1j, 1-1j, 1+1j))
    elif mod_size == 16:
        symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j, -1+1j,
                            -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j, 3-1j, 3+1j,
                            3+3j))
    signal_digmod = np.random.choice(symbols, size=(num_sym, 1), replace=True)

    return signal_digmod
