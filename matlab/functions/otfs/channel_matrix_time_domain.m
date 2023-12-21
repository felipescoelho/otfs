function [G] = channel_matrix_time_domain(gs, M, N, delaySpread)
% CHANNEL_MATRIX_TIME_DOMAIN    Generates the time domain matrix for the
%                               OTFS system.
%   [G] = CHANNEL_MATRIX_TIME_DOMAIN(gs, M, N, delaySpread) .
%
%   - Paramters:
%       . 

G = zeros(N*M, N*M);
for q = 0:N*M-1
    for ell = 0:delaySpread
        if q >= ell
            G(q+1, q+1-ell) = gs(ell+1, q+1);
        end
    end
end
end

