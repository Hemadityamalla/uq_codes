close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 1e2;

%f = @(x) 1.0*(x > 0.5);%Test function to be integrated
f = @(x) x.*(x < 0.5) + (x >= 0.5).*exp(-0.5*x);
%f = @(x) exp(-abs(x - 0.5)/2)/sqrt(2*pi);
%f = @(x) exp(-x.^2/2/sqrt(2*pi));
N = 5:5:50;
error = [];
for i = N
    %[xt,wt,samples] = almost_adaptive_quad(i,f,Kmax);
    [xt,wt,samples] = cappedaccuracy_quad([i+10,20],f,Kmax);
    error(end+1) = abs(dotprod(f(xt'),wt) - mean(f(samples)))
    
end
semilogy(N, error,'-o')