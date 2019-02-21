clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

f = @(x) exp(-0.5*sum(abs(x-0.5),2));
%f = @(x) cos(2*pi*0.3 + 0.5*(sum(x,2)));
%f = @(x) (prod(x,2) <= 0)*0 + (prod(x,2) > 0).*exp(0.5*sum(x,2)); %Discontinuous
d = 5; %dimension of the random vector
% Input for type of orthogonal polynomial basis. Alternatively, use 'Hermite'
polyBasis = 'Legendre';
quadrature = 'ClenshawCurtis';
growth = @(x) 2^(x-1);
exact = 0.5*(((1.0 - exp(-0.25))/0.5)^d + ((1.0 - exp(-0.25))/0.5)^d);%0.5*(((1.0 - exp(-0.25))/0.5)^d + ((1.0 - exp(-0.25))/0.5)^d);%cos(2*pi*0.3 + 0.25*d)*(sin(0.25)/0.5)^d;%((exp(0.5) - 1))^d;
error = [];
numpts = [];
for k=[1,2,3,4,5,6]

[eval_pts,weights] = smolyakSparseGrid(d,k,growth, quadrature);
%Rescaling to the domain [0,1].
eval_pts = 0.5*(eval_pts+1);
weights = weights/2^d;
%Quadrature rule for the mean square error(Q^d points)
%Q = 2;
%[eval_pts_mse,weights_mse] = smolyakSparseGrid(d, Q, growth, quadrature);

N = 0; %maximum degree of the multivariate polynomial
%MSE = []; %Empty array to store the mean-squared error
     
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
        %Integral = sum(legendre(eval_pts',lexOrdering(i_P,:)').*f(eval_pts').*prod(weights',2))
        fhat(i_P,1) = (sum(legendre(eval_pts',lexOrdering(i_P,:)').*f(eval_pts').*prod(weights',2)))/prod(gamma(lexOrdering(i_P,:)+1))
        %fapprox = fapprox + fhat(i_P,1)*legendre(eval_pts_mse', lexOrdering(i_P,:)');
    end
    error = [error;abs(fhat(1,1) - exact)];
    numpts = [numpts;length(weights)];
    %Computing the mean squared error
    %MSE = [MSE;sqrt(sum((fapprox - f(eval_pts_mse')).^2.*prod(weights_mse',2)))];
end

