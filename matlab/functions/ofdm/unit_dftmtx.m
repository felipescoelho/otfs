function [W] = unit_dftmtx(N)
% UNIT_DFTMTX   Generates the unitary DFT matrix.
%   [W] = UNIT_DFTMTX(N) .
%
%   - Parameters:
%       . N : Length of the DFT.
%   - Returns:
%       . W : (N, N) array with elements of the DFT matrix.
%

W = zeros(N, N);
for k = 1:N
    for l = 1:N
        W(k, l) = (1/sqrt(N))*exp(2*1j*pi*l*k/N);
    end
end
end
