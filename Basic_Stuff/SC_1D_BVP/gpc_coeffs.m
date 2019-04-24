function coeffs = gpc_coeffs(degree, nodes, weights, fn, type)
coeffs = zeros(degree+1,1);
if strcmp(type, 'hermite')
    gamma = factorial(0:degree);
    for ii=1:degree+1
       coeffs(ii,1) = sum(weights.*hermite(nodes, ii-1).*fn)/gamma(ii); 
    end
    
elseif strcmp(type, 'legendre')
    gamma = 2.0./(2*(0:degree) + 1.0);
    for ii=1:degree+1
       coeffs(ii,1) = dotprod(weights.*legendre(nodes, ii-1).*fn)/gamma(ii); 
    end
else
    error('Wrong basis!');
end
end