function [H, h] = channel_matrix_ofdm(gi, ki, li, numSubcarriers, numBlocks, cpLength)
% CHANNEL_MATRIX_OFDM   Generates the channel matrix for OFDM system with
%                       adjustable cyclic prefix.
%   [H, h] = CHANNEL_MATRIX_OFDM(gi, ki, li)
%
%   - Parameters:
%       . gi : Channel coefficients.
%       . li : Normalized delay taps.
%       . ki : Normalized Doppler taps.
%       . numSubcarriers : Length of the DFT.
%       . numBlocks : Number of information blocks.
%       . cpLength : Length for the cyclic prefix.
%   - Returns:
%       . H : Channel matrix.
%       . h : Channel TDL.
%

z = exp(1j*2*pi/(numSubcarriers+cpLength)*numBlocks);
delaySpread = max(li);

h = zeros(delaySpread+1, (numSubcarriers+cpLength)*numBlocks);
for q = 0:(numSubcarriers+cpLength)*numBlocks-1
    for i = 1:length(li)
        h(li(i)+1, q+1) = h(li(i)+1, q+1) + gi(i)*z^(ki(i) * (q-li(i)));
    end
end

H = zeros((numSubcarriers+cpLength)*numBlocks, (numSubcarriers+cpLength)*numBlocks);
for q = 0:(numSubcarriers+cpLength)*numBlocks-1
    for ell = 0:delaySpread
        if q >= ell
            H(q+1, q-ell+1) = h(ell+1, q+1);
        end
    end
end

end

