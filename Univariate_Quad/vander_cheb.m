function V = vander_cheb(N, pts)

%The below code can be vectorized
    V = zeros(N,N);
    for idx=1:N
        V(:,idx) = chebyshev(pts, idx-1);
    end
end