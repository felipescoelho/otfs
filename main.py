"""main.py

Script to run OFTS simulation.
Based on Delay-Doppler Communications: Principles and Applications -- Yi
Hong, Tharaj Thaj, and Emanuele Viterbo.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np
from src import dft_mat, qam_sym, qam_sym_decision


if __name__ == '__main__':
    # Definitions:
    # Frame parameters
    num_doppler_bins = 16  # N
    num_delay_bins = 64  # M
    dft_mtx = dft_mat(num_doppler_bins)
    dft_mtx = dft_mtx/np.linalg.norm(dft_mtx)  # Normalization
    delta_freq = 15*1e3
    block_duration = 1/delta_freq
    carrier_freq = 4*1e9
    delay_resol = 1/(num_delay_bins*delta_freq)
    doppler_resol = 1/(num_doppler_bins*block_duration)
    mod_size = 4
    
    # Transmitter:
    # OTFS frame
    num_sym_per_frame = num_delay_bins*num_doppler_bins
    tx_dig_sym = qam_sym(num_sym_per_frame, mod_size)
    # dely-Doppler frame
    x_frame = tx_dig_sym.reshape((num_delay_bins, num_doppler_bins))
    # Using IDZT:
    x_tilda_frame = x_frame @ np.conj(dft_mtx.T)  # Eq. 4.17
    s = x_tilda_frame.reshape((1, num_sym_per_frame))  # Eq. 4.20

    # Channel:
    

    # Receiver:
    r = s.copy()
    y_tilda_frame = r.reshape((num_delay_bins, num_doppler_bins))  # Eq. 4.25
    y_frame = y_tilda_frame @ dft_mtx  # Eq. 4.27

    rx_dig_sym = qam_sym_decision(y_frame.reshape((num_sym_per_frame, 1)),
                                  mod_size)
    
    print(sum(rx_dig_sym != tx_dig_sym))

    


