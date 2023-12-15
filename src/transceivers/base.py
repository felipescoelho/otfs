"""transceiver_base.py

Script with base transceiver.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


from abc import ABC, abstractclassmethod
from scipy.constants import speed_of_light
import numpy as np


class TranscieverBase(ABC):
    """"""

    @staticmethod
    def qam_sym(num_sym: int, mod_size: int):
        """Method to generate random QAM symbols.
        
        Parameters
        ----------
        num_sym : int
            Number of symbols.
        mod_size : int
            Size of the QAM constallation.
        
        Returns
        -------
        signal_digmod : np.ndarray
        """

        if mod_size == 4:
            symbols = np.array((-1-1j, -1+1j, 1-1j, 1+1j))
        elif mod_size == 16:
            symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j,
                                -1+1j, -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j,
                                3-1j, 3+1j, 3+3j))
        signal_digmod = np.random.choice(symbols, size=(num_sym, 1),
                                         replace=True)
        return signal_digmod
    
    @staticmethod
    def qam_sym_decision(x_hat: np.ndarray, mod_size: int):
        """Method to decide which symbol in rx.

        Parameters
        ----------
        x_hat : np.ndarray
            Received signal.
        mod_size : int
            Size of modulation.
        
        Returns
        -------
        rx_signal_digmod : np.ndarray
            Received information symbols.
        """
        if mod_size == 4:
            symbols = np.array((-1-1j, -1+1j, 1-1j, 1+1j))
        elif mod_size == 16:
            symbols = np.array((-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j,
                                -1+1j, -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j,
                                3-1j, 3+1j, 3+3j))
        num_sym = len(x_hat)
        rx_signal_digmod = np.zeros((num_sym, 1), dtype=np.complex128)
        for idx in range(num_sym):
            min_idx = np.argmin(np.abs(symbols - x_hat[idx]))
            rx_signal_digmod[idx] = symbols[min_idx]
        return rx_signal_digmod
    
    @abstractclassmethod
    def tx_to_rx(self):
        """Method for channel effects.
        
        Parameters
        ----------
        x : np.ndarray
            Transmitted signal, time domain
        channel_mat : np.ndarray
            Time domain channel matrix.
        snr_dB : float
            Signal to noise ratio in dB.
        
        Returns
        -------
        y : np.ndarray
            Received signal.
        sigma_n_2 : np.ndarray
            Noise power.
        """
    
    @abstractclassmethod
    def modulation(self):
        """Method for modulation.
        
        Parameters
        ----------
        x : np.ndarray
            Information symbols, in frame format MxN.

        Returns
        -------
        u : np.ndarray
            Modulation symbols to be transmitted, in column vector
            format MNx1.
        """

    @abstractclassmethod
    def demodulation(self):
        """Method for demodulation.
        
        Parameters
        ----------
        y : np.ndarray
            Received signal.
        
        Returns
        -------
        x_hat : np.ndarray
            Estimated symbols.
        """
    
    @abstractclassmethod
    def gen_channel(self):
        """"""

    def gen_channel_params(self, max_speed: float, doppler_resol: float,
                           delay_resol: float, model: str):
        """"""
        max_speed_ms = max_speed * (1000/3600)
        max_doppler = (max_speed_ms*self.carrier_freq)/speed_of_light
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
            g_i = np.sqrt(pdp)*(np.sqrt(.5)*(np.random.randn(num_taps,)
                                            + 1j*np.random.randn(num_taps,)))
            delay_taps = np.round(delays/delay_resol)
            doppler_taps = max_doppler_normalized*np.cos(
                2*np.pi*np.random.randn(num_taps,)
            )  # Assumes Jakes Spectrum
        elif model == 'synthetic':
            num_taps = 6
            max_delay_normalized = 4
            max_doppler_normalized = 4
            g_i = np.sqrt(1/num_taps)*(
                np.sqrt(.5) * (np.random.randn(num_taps,)
                            + 1j*np.random.randn(num_taps,))
            )
            delay_taps = np.random.randint(0, max_delay_normalized, (1, num_taps))
            delay_taps -= np.min(delay_taps)
            doppler_taps = max_doppler_normalized \
                - 2*max_doppler_normalized*np.random.rand(num_taps,)
        return g_i, delay_taps, doppler_taps

    def __init__(self, M: int, N: int, delta_freq: float, carrier_freq: float,
                 mod_size:int):
        """"""
        
        self.M = M
        self.N = N
        self.delta_freq = delta_freq
        self.carrier_freq = carrier_freq
        self.mod_size = mod_size

    def run_ser(self, mc: int, snr: np.ndarray, chann: np.ndarray):
        """Method to run simulation using Monte Carlo."""

        len_snr = len(snr)
        ser = np.zeros((len_snr,), dtype=np.float64)
        for idx in range(len_snr):
            for _ in range(mc):
                tx_dig_sym = self.qam_sym(int(self.M*self.N), self.mod_size)
                x_frame = tx_dig_sym.reshape(self.M, self.N)
                u = self.modulation(x_frame)
                y = self.tx_to_rx(u, chann, snr[idx])
                x_hat = self.demodulation(y)
                rx_dig_sym = self.qam_sym_decision(x_hat.flatten(), self.mod_size)
                ser[idx] += sum(rx_dig_sym != tx_dig_sym)
            ser[idx] /= mc
        return ser
