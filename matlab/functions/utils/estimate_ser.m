function [ser] = estimate_ser(xHat, txSymbols, modSize)
% ESTIMATE_SER  Function to estimate symbol error rate.
%   [ser] = ESTIMATE_SER(xHat, txSymbols, modSize)
%
%   - Parameters:
%       . xHat : An estimate of the received symbols.
%       . txSymbols : Transmitted symbols.
%       . modSize : QAM constellation size.
%   - Returns:
%       . ser : Estimated symbol erro rate.
%

switch modSize
    case 4
        symbols = [-1-1j, -1+1j, 1-1j, 1+1j];
    case 16
        symbols = [-3-3j, -3-1j, -3+1j, -3+3j, -1-3j, -1-1j, -1+1j, ...
            -1+3j, 1-3j, 1-1j, 1+1j, 1+3j, 3-3j, 3-1j, 3+1j, 3+3j];
end

xHatVect = xHat(:);
rxSymbolsEstimated = zeros(length(xHatVect), 1);
for idx = 1:length(xHatVect)
    [~, symbolIdx] = min(sqrt(symbols - xHatVect(idx)));
    rxSymbolsEstimated(idx) = symbols(symbolIdx);
end
ser = sum(rxSymbolsEstimated ~= txSymbols(:))/length(xHatVect);
end
