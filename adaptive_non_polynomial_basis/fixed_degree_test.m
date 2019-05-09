close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 1e3;

%f = @(x) 1.0*(x > 0.5);%Test function to be integrated
%f = @(x) x.*(x < 0.5) + (x >= 0.5).*exp(-0.5*x);
f = @(x) exp(-abs(x - 0.5)/2)/sqrt(2*pi);
%f = @(x) exp(-x.^2/2/sqrt(2*pi));
N = 5:5:45;
error = [];
for i = N
    [xt,wt,samples] = fixed_implict_quad(i,Kmax);
    error(end+1) = abs(dotprod(f(xt'),wt) - mean(f(samples)))
    
end
semilogy(N, error,'-o');