function [R] = rm_redundancy_matrix(dftLength, cpLength)
% RM_REDUNDANCY_MATRIX   Generates matrix to remove redundancy from CP-OFDM
%                       frame.
%   [R] = RM_REDUNDANCY_MATRIX(dftLength, cpLength)
%
%   - Parameters:
%       . dftLength : Length of the DFT (number of subcarriers).
%       . cpLength : Number of samples in cyclic prefix.
%   - Returns:
%       . R : Matrix that removes redundancy.
%

R = [zeros(dftLength, cpLength), eye(dftLength)];
end

