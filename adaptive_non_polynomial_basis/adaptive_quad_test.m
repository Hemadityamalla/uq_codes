close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 1e4;

f = @(x) (x > 0.5).*(exp(0.5*x));%Test function to be integrated
%f = @(x) exp(-x.^2/2)/sqrt(2*pi);
N = 2:4:30;
error = [];
for i = N
    [xt,wt,samples] = almost_adaptive_quad(i,f,Kmax);
    error(end+1) = abs(dotprod(f(xt'),wt) - mean(f(samples)))
    
end
semilogy(N, error,'-o')