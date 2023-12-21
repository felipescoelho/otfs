function [y, noisePower] = add_wgn(x, snr)
% ADD_WGN   Adds white Gaussian noise with adjuted power by SNR in dB.
%   [y, noisePower] = ADD_WGN(x, snr) .
%
%   - Parameters:
%       . x : Signal to be corrupted with noise.
%       . snr : Noise power in dB.
%   - Return:
%       . y : Noisy signal.
%       . noisePower : Noise power according to SNR and signal noise.
%

if isreal(x)
    noisePre = randn(size(x))+1j*randn(size(x));
else
    noisePre = randn(size(x));
end
noisePowerPre = mean(abs(noisePre).^2);
noisePower = 10^(-.1*snr)*mean(abs(x).^2)/noisePowerPre;
noise = sqrt(noisePower)*noisePre;
y = x+noise;
end
