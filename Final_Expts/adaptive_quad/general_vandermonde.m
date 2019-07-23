function V = general_vandermonde(pts, basis, params)


%The below code can be vectorized
    V = zeros(length(params), length(pts));
    for idx=1:length(params)
        V(idx,:) = basis(pts, params(idx));
    end
    end