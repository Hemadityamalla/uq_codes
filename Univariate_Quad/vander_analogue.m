function V = vander_analogue(N, pts, fn)

%The below code can be vectorized
    V = zeros(N,N);
    for idx=1:N
        V(:,idx) = fn(pts, idx-1);
    end
end