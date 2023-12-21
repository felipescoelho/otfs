% otfs_ber.m
%
% OTFS -- BER simulation for different values of SNR using Monte Carlo.
%
% Dec. 17, 2023
%

clc
clear
close all


addpath('functions')


% OTFS - Settings:
N = 16;  % Number of Doppler bins
M = 64;  % Number of delay bins
Fn = dftmtx(N);
Fn = Fn/norm(Fn);  % Normalized DFT
deltaF = 15*1e3;  % Subcarrier spacing
T = 1/deltaF;  % Block duration
fc = 4*1e9;  % Carrier frequency
c = physconst('Lightspeed');
delayResolution = 1/(M*deltaF);
dopplerResolution = 1/(N*T);
modSize = 4;
numSymbolsPerFrame = N*M;
numBitsPerFrame = numSymbolsPerFrame*log2(modSize);

% Channel EVA
maxUserSpeed = 250*(1000/3600);  % m/s
nuMax = maxUserSpeed*fc/c;
kMax = nuMax/dopplerResolution;
delays = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510]*1e-9;
pdp = [0, -1.5, -1.4, -3.6, -.6, -9.1, -7, -12, -16.9];
pdpLinear = 10.^(pdp/10);
pdpLinear = pdpLinear/sum(pdpLinear);
numTaps = length(pdpLinear);
gi = sqrt(pdpLinear).*(sqrt(.5)*(randn(1, numTaps) + 1i*randn(1, numTaps)));
li = round(delays./delayResolution);  % Delay taps
ki = (kMax*cos(2*pi*randn(1, numTaps)));
z = exp(1i*2*pi/(N*M));
delaySpread = max(li);
gs = zeros(delaySpread+1, N*M);
for q = 0:N*M-1
    for i = 1:numTaps
        gs(li(i)+1, q+1) = gs(li(i)+1, q+1) + gi(i)*z^(ki(i)*(q-li(i)));
    end
end
G = zeros(N*M, N*M);
for q = 0:N*M-1
    for ell = 0:delaySpread
        if q >= ell
            G(q+1, q-ell+1) = gs(ell+1, q+1);
        end
    end
end

% Simulation settings:
monteCarloEnsemble = 50;
snrdB = linspace(-15, 40, 25);

% Memory allocation:
BERNoise = zeros(1, length(snrdB));
BERTimeDomainLMMSE = zeros(1, length(snrdB));
BERGuilherme = zeros(1, length(snrdB));
meanBERGuilherme = zeros(1, length(snrdB));

for snrIdx = 1:length(snrdB)
    ber = 0;
    berNoise = 0;
    for it = 1:monteCarloEnsemble
        txBits = randi([0, 1], numBitsPerFrame, 1);
        txSymbols = qammod(txBits, modSize, 'gray', 'InputType', 'bit');
        X = reshape(txSymbols, M, N);
        XTilda = X*Fn';
        s = reshape(XTilda, 1, N*M);
        % Time domain channel matrix
        r = G*s.';
        % Add thermal noise, and get power estimate:
        snr = 10.^(snrdB(snrIdx)/10);
        Es = mean(abs(r).^2);
        sigmaN2 = Es/snr;
        noise = sqrt(sigmaN2/Es)*(randn(N*M, 1) + 1j*randn(N*M, 1));
        r = r+noise;
        
        % Time domain LMMSE detection
        sHat = (G'*G + sigmaN2*eye(M*N)) \ (G'*r);
        XHatTilda = reshape(sHat, M, N);
        XHat = XHatTilda*Fn;
        xHat = reshape(XHat, N*M, 1);
        rxBits = qamdemod(xHat, modSize, 'gray', 'OutputType', 'bit');
        err = sum(rxBits ~= txBits);
        ber = ber + err/length(txBits);
        BERGuilherme(1, snrIdx) = BERGuilherme(1, snrIdx) + err;
        meanBERGuilherme(1, snrIdx) = BERGuilherme(1, snrIdx)*1/length(txBits)*1/it;
        
        % AWGN channel
        rNoise = s.'+noise;
        XHatTildaNoise = reshape(rNoise, M, N);
        XHatNoise = XHatTildaNoise*Fn;
        xHatNoise = reshape(XHatNoise, N*M, 1);
        rxBitsNoise = qamdemod(xHatNoise, modSize, 'gray', ...
            'OutputType', 'bit');
        errNoise = sum(rxBitsNoise ~= txBits);
        berNoise = berNoise + errNoise/length(txBits);
    end
    BERTimeDomainLMMSE(1, snrIdx) = ber/monteCarloEnsemble;
    BERNoise(1, snrIdx) = berNoise/monteCarloEnsemble;
end

figName = strjoin([{'BER OTFS -- Max speed: '}, ...
    num2str(maxUserSpeed*3.6), {' km/h.'}]);
fig = figure;
fig.Name = figName;
semilogy(snrdB, BERTimeDomainLMMSE, '-s'), hold on
semilogy(snrdB, meanBERGuilherme, '-*')
semilogy(snrdB, BERNoise, '-d'), hold off, grid on
legend('BER - Time Domain Matrix LMMSE', 'BER - Only Noise')
ylim([1e-6, 1])
xlabel('SNR, dB')
ylabel('BER')


% EoF
