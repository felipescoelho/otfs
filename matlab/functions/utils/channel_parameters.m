function [gi, ki, li] = channel_parameters(channelMode, maxSpeed, carrierFreq, delayResolution, dopplerResolution)
% CHANNEL_PARAMETERS    Generates set of parameters to simulate multipath
%                       channels for OFDM and OTFS.
%   [gi, ki, li] = CHANNEL_PARAMETERS(channelMode, maxSpeed, carrierFreq,
%                                     delayResolution, dopplerResolution) .
%
%   - Paramters:
%       . channelMode : Which type of channel is used.
%       . maxSpeed : Maximum relative speed, in km/h, between Tx and Rx.
%       . carrierFreq : Carrier frequency.
%       . delayResolution : Sampling period for delay axis. The same as the
%                           transmitted signal sampling period.
%       . dopplerResolution : Sampling period for Doppler axis. The same as
%                             subcarrier carrier spacing divided by the
%                             number of transmitted blocks.
%   - Returns:
%       . gi : Channel coefficients.
%       . li : Delay taps.
%       . ki : Doppler taps.
%

possibleMode = {'EPA', 'EVA', 'ETU'};
errorMessage = 'The possible channel models are EPA, EVA, or ETU.';
assert(ismember(channelMode, possibleMode), errorMessage)

maxSpeedmps = maxSpeed/3.6;
maxDoppler = maxSpeedmps*carrierFreq/physconst('Lightspeed');
switch channelMode
    case 'EPA'
        delays = [0, 30, 70, 90, 110, 190, 410] * 1e-9;
        pdp = [0, -1, -2, -3, -8, -17.2, -20.8];
    case 'EVA'
        delays = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
        pdp = [0, -1.5, -1.4, -3.6, -.6, -9.1, -7, -12, -16.9];
    case 'ETU'
        delays = [0, 50, 120, 200, 230, 500, 1600, 2300, 5000] * 1e-9;
        pdp = [-1, -1, -1, 0, 0, 0, -3, -5, -7];
end
pdpLinear = 10.^(.1*pdp);
pdpLinear = pdpLinear/sum(pdpLinear);
numPaths = length(pdpLinear);
complexNoise = randn(numPaths, 1) + 1j*randn(numPaths, 1);
gi = sqrt(pdpLinear/mean(abs(complexNoise).^2)).*complexNoise;
li = round(delays/delayResolution);
ki = maxDoppler*cos(2*pi*rand(numPaths, 1))/dopplerResolution;
end

