clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

f = @(x) exp(sum(abs(x),2));
%f = @(x) cos(2*pi*0.3 + 0.5*(sum(x,2)));
%f = @(x) (x(:,1) > 0)*0 .* (x(:,2) > 0)*0 + (x(:,1) <= 0).*exp(0.5*(x(:,1) + x(:,2))).*(x(:,2) <= 0);
d = 15; %dimension of the random vector
% Input for type of orthogonal polynomial basis. Alternatively, use 'Hermite'
polyBasis = 'Legendre';
quadrature = 'ClenshawCurtis';
growth = @(x) 2^(x-1);
k = 3;
[eval_pts,weights] = smolyakSparseGrid(d,k,growth, quadrature);

%Quadrature rule for the mean square error(Q^d points)
Q = 2;
[eval_pts_mse,weights_mse] = smolyakSparseGrid(d, Q, growth, quadrature);

order = [0,1]; %maximum degree of the multivariate polynomial
MSE = []; %Empty array to store the mean-squared error
for N = order
     
    lexOrdering = monomialDegrees(d,N);
    %Pre-computing the normalization factors
    if strcmp(polyBasis,'Hermite')
        gamma = factorial(0:N); %Hermite
    else
        gamma = 2.0./(2*(0:N) + 1.0); %Legendre
    end
    P = size(lexOrdering,1);
    fhat = zeros(P,1);
    fapprox = 0;
    for i_P = 1:P
        fhat(i_P,1) = (sum(legendre(eval_pts',lexOrdering(i_P,:)').*f(eval_pts').*prod(weights',2)))/prod(gamma(lexOrdering(i_P,:)+1));
        fapprox = fapprox + fhat(i_P,1)*legendre(eval_pts_mse', lexOrdering(i_P,:)');
    end
    %Computing the mean squared error
    MSE = [MSE;sqrt(sum((fapprox - f(eval_pts_mse')).^2.*prod(weights_mse',2)))];
end

