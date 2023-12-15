"""channel.py

Script with method to generate channel for OTFS.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np
from scipy.constants import speed_of_light


def gen_chan_params(max_speed:float, carrier_freq:float, doppler_resol:float,
                    delay_resol:float, model:str):
    """Method to generate channel model parameters.
    
    Parameters
    ----------
    max_speed : float
        Maximum user equipment speed (km/h).
    carreier_freq : float
        Carrier frequency.
    doppler_resol : float
        Resolution for Doppler axis.
    model : str
        Selcted channel model.
    
    Returns
    -------
    g_i : np.ndarray
    delay_taps : np.ndarray
    doppler_taps : np.ndarray
    """

    max_speed_ms = max_speed * (1000/3600)
    max_doppler = (max_speed_ms*carrier_freq)/speed_of_light
    max_doppler_normalized = max_doppler/doppler_resol

    if model in ['EPA', 'EVA', 'ETU']:
        if model == 'EPA':
            delays = np.array((0, 30, 70, 90, 110, 190, 410))*1e-9
            pdp_db = np.array((0, -1, -2, -3, -8, -17.2, -20.8))
        elif model == 'EVA':
            delays = np.array((0, 30, 150, 310, 370, 710, 1090, 1730, 2510))*1e-9
            pdp_db = np.array((0, -1.5, -1.4, -3.6, -.6, -9.1, -7, -12, -16.9))
        elif model == 'ETU':
            delays = np.array((0, 50, 120, 200, 230, 500, 1600, 2300, 5000))*1e-9
            pdp_db = np.array((-1, -1, -1, 0, 0, 0, -3, -5, -7))
        pdp = 10**(pdp_db/10)
        num_taps = len(pdp)
        g_i = np.sqrt(pdp)*(np.sqrt(.5)*(np.random.randn(1, num_taps)
                                         + 1j*np.random.randn(1, num_taps)))
        delay_taps = np.round(delays/delay_resol)
        doppler_taps = max_doppler_normalized*np.cos(
            2*np.pi*np.random.randn(num_taps,)
        )  # Assumes Jakes Spectrum
    elif model == 'synthetic':
        num_taps = 6
        max_delay_normalized = 4
        max_doppler_normalized = 4
        g_i = np.sqrt(1/num_taps)*(
            np.sqrt(.5) * (np.random.randn(1, num_taps)
                           + 1j*np.random.randn(1, num_taps))
        )
        delay_taps = np.random.randint(0, max_delay_normalized, (1, num_taps))
        delay_taps -= np.min(delay_taps)
        doppler_taps = max_doppler_normalized \
            - 2*max_doppler_normalized*np.random.rand(1, num_taps)
    
    return g_i, delay_taps, doppler_taps


def gen_delay_time(num_sym_per_frame:int, delay_taps:np.ndarray, g_i:np.ndarray,
                   doppler_taps:np.ndarray):
    """Method to generate delay-time channel coefficientes and matrix.
    
    Parameters
    ----------
    num_sym_per_frame : int
        Number of symbols per frame.
    delay_taps : np.ndarray
        Delay taps for the channel.
    g_i : np.ndarray
        Channel coefficients with uniform pdp.
    doppler_taps: np.ndarray
        Doppler taps for the channel.
    
    Returns
    -------
    g_s : np.ndarray
        Channel in the for of TDL.
    G : np.ndarray
        Baseband channel matrix.
    delay_spread : float
        Delay spread.
    """

    z = np.exp(1j*2*np.pi/num_sym_per_frame)
    delay_spread = int(np.max(delay_taps))
    num_taps = len(delay_taps)
    # Channel TDL:
    gs = np.zeros((delay_spread+1, num_sym_per_frame), dtype=np.complex128)
    for q in range(num_sym_per_frame):
        for i in range(num_taps):
            delay_tap = int(delay_taps[i])
            gs[delay_tap, q] = gs[delay_tap, q] + g_i[0, i] \
                * z**(doppler_taps[i]*(q-delay_tap))
    # Baseband channel matrix:
    G = np.zeros((num_sym_per_frame, num_sym_per_frame), dtype=np.complex128)
    for q in range(num_sym_per_frame):
        for ell in range(delay_spread+1):
            if q >= ell:
                G[q, q-ell] = gs[ell, q]
    
    return gs, G, delay_spread


def gen_delay_time_and_delay_doppler(num_delay_bins:int, P_mat:np.ndarray,
                                     G:np.ndarray, dft_mtx:np.ndarray):
    """Method to generate delay-time and delay-Doppler channel matrix.
    
    Parameters
    ----------
    num_doppler_bins : int
        Number of Doppler bins.
    P_mat : np.ndarray
        Row-column permutation matrix.
    G : np.ndarray
        Discrete-time baseband channel matrix.
    dft_mtx : np.ndarray
        Normalized DFT matrix.

    Returns
    -------
    H_tilda : np.ndarray
        Delay-time channel matrix.
    H : np.ndarray
        Delay-Doppler channel matrix.
    """

    Im = np.eye(num_delay_bins)
    H_tilda = P_mat.T @ G @ P_mat
    H = np.kron(Im, dft_mtx) @ (np.conj(P_mat.T)@G@P_mat) \
        @ np.kron(Im, np.conj(dft_mtx.T))

    return H_tilda, H


def filter_sig(s:np.ndarray, method:str, *args):
    """Method to filter signal s.
    Generates r by passing the tx through the channel.
    
    Parameters
    ----------
    s : np.ndarray
        Input signal.
    method : str
        Method to filter signal. (m1, m2, m3, m4)
    *args : tuple
        Arguments according to selected method.
        if m1: gs, num_sym_per_frame, delay_spread
        elif m2: G
        elif m3: H_tilda, num_sym_per_frame, P_mat
        elif m4: H, num_sym_per_frame, P_mat, num_doppler_bins, dft_mtx
    
    Returns
    -------
    r : np.ndarray
        Received signal without noise.
    """

    if method == 'm1':
        gs, num_sym_per_frame, delay_spread = args
        r = np.zeros((num_sym_per_frame, 1), dtype=np.complex128)
        for q in range(num_sym_per_frame):
            for ell in range(delay_spread):
                if q >= ell:
                    r[q] += gs[ell, q]*s[q-ell]
    elif method == 'm2':
        G = args
        r = G @ s.T
    elif method == 'm3':
        H_tilda, num_sym_per_frame, P_mat = args
        x_tilda = (s.T).reshape((num_sym_per_frame, 1))
        y_tilda = H_tilda @ x_tilda
        r = P_mat @ y_tilda
    elif method == 'm4':
        H, num_sym_per_frame, P_mat, num_doppler_bins, dft_mtx = args
        Im = np.eye(num_doppler_bins)
        x = s.reshape((num_sym_per_frame, 1))
        y = H @ x
        r = P_mat @ np.kron(Im, np.conj(dft_mtx.T)) @ y 
        
    return r