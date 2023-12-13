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
            2*np.pi*np.random.randn(1, num_taps)
        )  # Assumes Jakes Spectrum
    elif model == 'synthetic':
        num_taps = 6
        max_delay_nomalized = 4
        max_doppler_normalized = 4
        g_i = np.sqrt(1/num_taps)*(
            np.sqrt(.5) * (np.random.randn(1, num_taps)
                           + 1j*np.random.randn(1, num_taps))
        )
        delay_taps = np.random.randint(0, max_delay_nomalized, (1, num_taps))
        delay_taps -= np.min(delay_taps)
        doppler_taps = max_doppler_normalized \
            - 2*max_doppler_normalized*np.random.rand(1, num_taps)
    
    return g_i, delay_taps, doppler_taps


def gen_delay_time(num_sym_per_frame:int, delay_taps:np.ndarray, g_i:np.ndarray,
                   doppler_taps:np.ndarray):
    """Method to generate delay-time channel coefficientes and matrix.
    
    Parameters
    ----------
    
    Returns
    -------
    """

    z = np.exp(1j*2*np.pi/num_sym_per_frame)
    delay_spread = np.max(delay_taps) - np.min(delay_taps)
    num_taps = len(delay_taps)
    # Channel TDL:
    gs = np.zeros((delay_spread+1, num_sym_per_frame), dtype=np.complex128)
    for q in range(num_sym_per_frame):
        for i in range(num_taps):
            gs[delay_taps[i], q] += g_i[i]*z**(doppler_taps[i]*1-delay_taps[i])
    # Baseband channel matrix:
    G = np.zeros((num_sym_per_frame, num_sym_per_frame), dtype=np.complex128)
    for q in range(num_sym_per_frame):
        for ell in range(delay_spread+1):
            if q >= ell:
                G[q, q-ell] = gs[ell, q]
    
