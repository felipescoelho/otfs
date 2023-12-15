"""ofdm.py

Script with class for OFDM transceivers.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np
from .base import TranscieverBase


class OFDM(TranscieverBase):
    """"""

    def __init__(self, no_symbols: int, dft_len: int, delta_freq: float,
                 carrier_freq: float, mod_size: int, cp_len: int):
        super().__init__(no_symbols, dft_len, delta_freq, carrier_freq, mod_size)
        self.cp_len = cp_len
    
    def add_redundancy_matrix(self):
        return np.vstack((
            np.hstack((np.zeros((self.cp_len, self.N-self.cp_len)),
                       np.eye(self.cp_len))),
            np.eye(self.N)
        ))
    
    def rm_redudancy_matrix(self):
        return np.hstack((np.zeros((self.N, self.cp_len)), np.eye(self.N)))
    
    def dft_matrix(self):
        dft_mat = np.zeros((self.N, self.N), dtype=np.complex128)
        for k in range(self.N):
            for l in range(self.N):
                dft_mat[k, l] = np.exp(-1j*2*np.pi*k*l/self.N)
        return dft_mat
    
    def idft_matrix(self):
        idft_mat = np.zeros((self.N, self.N), dtype=np.complex128)
        for k in range(self.N):
            for l in range(self.N):
                idft_mat[k, l] = np.exp(1j*2*np.pi*k*l/self.N)
        idft_mat /= self.N
        return idft_mat

    def modulation(self, x: np.ndarray):
        return (self.add_redundancy_matrix()@self.idft_matrix()@ x.T)\
            .reshape((int(self.M*(self.N+self.cp_len)), 1))
    
    def demodulation(self, y: np.ndarray):
        return (self.dft_matrix()@self.rm_redudancy_matrix()
                @ y.reshape((self.N+self.cp_len, self.M))).T
    
    def gen_channel(self, max_speed: float, doppler_resol: float,
                    delay_resol: float):
        g_i, delay_taps, doppler_taps = self.gen_channel_params(
            max_speed, doppler_resol, delay_resol
        )
        z = np.exp(1j*2*np.pi/(self.N*self.M))
        delay_spread = int(max(delay_taps))
        gs = np.zeros((delay_spread+1, self.N*(self.M+self.cp_len)))
        for q in range(self.N*(self.M+self.cp_len)):
            for i in range(len(delay_taps)):
                gs[int(delay_taps[i]), q] += g_i[i] \
                    * z**(doppler_taps[i]*(q-int(delay_taps[i])))
                
        
    def tx_to_rx(self, x: np.ndarray, channel_mat: np.ndarray, snr_dB: float):
        r = channel_mat @ x
        Es = np.mean(np.abs(r)**2)
        snr = 10**(snr_dB/10)
        sigma_n_2 = Es/snr
        noise = np.sqrt(sigma_n_2/2)*(np.random.randn(len(r), 1)
                                      +1j*np.random.randn(len(r), 1))
        r += noise
        return r, sigma_n_2

