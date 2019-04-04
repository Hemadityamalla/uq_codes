function v = haar(N, x)
% HAAR(N, x) - Evaluate N-th Haar basis function on x
%
% Implementation is really inefficient.

    i = 0;
    for n=0:N
        for k=0:2^n-1
            if i == N
                v = psi(x, n, k);
                return;
            end
            i = i+1;
        end
    end
end

% Wavelet mother function
function v = psi(t, n, k)
    if nargin < 2
        v = t;
        v((0 <= t) & (t < 1/2)) = 1;
        v((1/2 <= t) & (t < 1)) = -1;
        v((t < 0) | (t >= 1)) = 0;
    else
        v = psi(2^n*t-k);
    end
end