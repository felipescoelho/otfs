function [x] = qam_random_symbols(modSize, M, N)
% QAM_RANDOM_SYMBOLS    Generates random symbols from QAM constellation
%                       with modulation size defined by modSize.
%   [x] = QAM_RANDOM_SYMBOLS(modSize, M, N) .
%
%   - Parameters:
%       . modSize : QAM constellation size.
%       . M : Number of rows in x.
%       . N : Number of columns in x.
%   - Returns:
%       . x : (M, N) array with QAM symbols.
%

possibleModSize = [4, 16];
errorMessage = ['Select modSize values as' repmat(' %d,', 1, ...
    numel(possibleModSize)-1) ' or %d.'];
assert(ismember(modSize, possibleModSize), errorMessage, ...
    possibleModSize(1), possibleModSize(2))

switch modSize
    case 4
        symbols = [-1-1j, -1+1j, 1-1j, 1+1j];
    case 16
        symbols = [-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j, -1+1j, ...
            -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j, 3-1j, 3+1j, 3+3j];
end

x = randsample(symbols, M*N, true);
x = reshape(x, M, N);
end
