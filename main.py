"""main.py

Script to run OFTS simulation.
Based on Delay-Doppler Communications: Principles and Applications -- Yi
Hong, Tharaj Thaj, and Emanuele Viterbo.

luizfelipe.coelho@smt.ufrj.br
Dec 13, 2023
"""


import numpy as np
from src import (dft_mat, qam_sym, qam_sym_decision,
                 gen_delay_time_and_delay_doppler, gen_delay_time,
                 gen_chan_params, filter_sig, row_col_permutation_mat)
# from src import OFDM, gen_chan_params, gen_delay_time

if __name__ == '__main__':

#     # Definitions
#     M = 64
#     N = 64
#     delta_freq = 15*1e3
#     block_duration = 1/delta_freq
#     carrier_freq = 4*1e9
#     delay_resol = 1/(M*delta_freq)
#     doppler_resol = 1/(N*block_duration)
#     mod_size = 4
#     snr_arr = np.linspace(-10, 50, 25)
#     # Models
#     ofdm_model = OFDM(no_symbols=M, dft_len=N, delta_freq=delta_freq,
#                       carrier_freq=carrier_freq, mod_size=mod_size, cp_len=16)
#     # Generate channels
#     max_speed_kmh = 320
#     g_i, delay_taps, doppler_taps = gen_chan_params(
#         max_speed_kmh, carrier_freq, doppler_resol, delay_resol, model='EVA'
#     )
#     _, G, _ = gen_delay_time(int(N*M), delay_taps, g_i,
#                                           doppler_taps)
#     # Run simulation
#     print(ofdm_model.run_ser(1, snr_arr, G))

    # Definitions:
    # Frame parameters
    num_delay_bins = 64  # M
    num_doppler_bins = 16  # N
    dft_mtx = dft_mat(num_doppler_bins)
    dft_mtx = dft_mtx/np.linalg.norm(dft_mtx)  # Normalization
    delta_freq = 15*1e3
    block_duration = 1/delta_freq
    carrier_freq = 4*1e9
    delay_resol = 1/(num_delay_bins*delta_freq)
    doppler_resol = 1/(num_doppler_bins*block_duration)
    mod_size = 4
    snr_db = 10
    
    # Transmitter:
    P_mat = row_col_permutation_mat(num_doppler_bins, num_delay_bins)
    # OTFS frame
    num_sym_per_frame = int(num_delay_bins*num_doppler_bins)
    tx_dig_sym = qam_sym(num_sym_per_frame, mod_size)
    # Dely-Doppler frame
    x_frame = tx_dig_sym.reshape((num_delay_bins, num_doppler_bins))
    # Using IDZT:
    x_tilda_frame = x_frame @ np.conj(dft_mtx.T)  # Eq. 4.17
    s = x_tilda_frame.reshape((1, num_sym_per_frame))  # Eq. 4.20

    # Channel:
    max_speed_kmh = 320
    g_i, delay_taps, doppler_taps = gen_chan_params(
        max_speed_kmh, carrier_freq, doppler_resol, delay_resol, model='EVA'
    )
    g_s, G, delay_spread = gen_delay_time(num_sym_per_frame, delay_taps, g_i,
                                          doppler_taps)
    H_tilda, H = gen_delay_time_and_delay_doppler(num_delay_bins, P_mat, G,
                                                  dft_mtx)
    args = (G)
    r = filter_sig(s, 'm2', *args)

    # Receiver:
    Es = np.mean(np.abs(tx_dig_sym)**2)
    snr = 10**(snr_db/10)
    sigma_w_2 = Es/snr
    noise = np.sqrt(sigma_w_2/2)*(np.random.randn(num_sym_per_frame, 1)
                                  + 1j*np.random.randn(num_sym_per_frame, 1))
    r += noise
    # OTFS demod:
    y_tilda_frame = r.reshape((num_delay_bins, num_doppler_bins))  # Eq. 4.25
    y_frame = y_tilda_frame @ dft_mtx  # Eq. 4.27
    rx_dig_sym_no_lmmse = qam_sym_decision(y_frame.reshape((num_sym_per_frame, 1)),
                                           mod_size)
    # Delay-Doppler LMMSE detection:
    y_dd = y_frame.reshape((num_sym_per_frame, 1))
    x_hat_dd = np.linalg.inv(np.conj(H.T)@H + sigma_w_2) @ (np.conj(H.T)@y_dd)
    rx_dig_sym_dd_lmmse = qam_sym_decision(x_hat_dd, mod_size)
    # Time domain LMMSE detection:
    s_hat = np.linalg.inv(np.conj(G.T)@G + sigma_w_2) @ (np.conj(G.T)@r)
    x_hat_tilda_frame = s_hat.reshape((num_delay_bins, num_doppler_bins))
    x_hat_frame = x_hat_tilda_frame @ dft_mtx
    x_hat_td = x_hat_frame.reshape((num_sym_per_frame, 1))
    rx_dig_sym_td_lmmse = qam_sym_decision(x_hat_td, mod_size)
    
    print(sum(rx_dig_sym_no_lmmse != tx_dig_sym))
    print(sum(rx_dig_sym_dd_lmmse != tx_dig_sym))
    print(sum(rx_dig_sym_td_lmmse != tx_dig_sym))
