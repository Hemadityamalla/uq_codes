function coeffs = gpc_coeffs(degree, nodes, weights, fn, type)
if size(degree,2) == 1 %Univariate case
    coeffs = zeros(degree+1,1);
    if strcmp(type, 'hermite')
        gamma = factorial(0:degree);
        for ii=0:degree
           coeffs(ii,1) = sum(weights.*hermite(nodes, ii).*fn)/gamma(ii); 
        end

    elseif strcmp(type, 'legendre')
        gamma = 2.0./(2*(0:degree) + 1.0);
        for ii=0:degree
           coeffs(ii,1) = sum(weights.*legendre(nodes, ii).*fn)/gamma(ii); 
        end
    end
else %Multivariate case
    P = size(degree,1);
    maxDegree = max(max(degree));
    coeffs = zeros(P,1);
    if strcmp(type, 'hermite')
        gamma = factorial(0:maxDegree);
        for ii=1:P
           coeffs(ii,1) = sum(weights.*hermite(nodes, degree(ii,:)').*fn)/prod(gamma(degree(ii,:)+1)); 
        end

    elseif strcmp(type, 'legendre')
        gamma = 2.0./(2*(0:maxDegree) + 1.0);
        for ii=1:P
           coeffs(ii,1) = sum(weights.*legendre(nodes, degree(ii,:)').*fn)/prod(gamma(degree(ii,:)+1)); 
        end
    end
    
end
    
end