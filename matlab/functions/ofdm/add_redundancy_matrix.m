function [Gamma] = add_redundancy_matrix(dftLength, cpLength)
% ADD_REDUNDANCY_MATRIX     Generates matrix to add redundancy in the form
%                           of cyclic prefix.
%   [Gamma] = ADD_REDUNDANCY_MATRIX(dftLength, cpLength)
%
%   - Parameters:
%       . dftLength : Length of the DFT (number of subcarriers).
%       . cpLength : Number of samples in cyclic prefix.
%   - Returns:
%       . Gamma : Matrix that adds redundancy to each block in a frame.
%

Gamma = [zeros(cpLength, dftLength-cpLength), eye(cpLength); eye(dftLength)];
end

