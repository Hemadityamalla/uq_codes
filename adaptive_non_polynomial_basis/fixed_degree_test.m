close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 1e4;

f = @(x) (x > 0.5).*(exp(0.5*x));%Test function to be integrated
N = 2:4:14;
error = [];
for i = N
    [xt,wt,samples] = fixed_implict_quad(i,Kmax);
    error(end+1) = abs(dotprod(f(xt'),wt) - 2*(exp(0.5) - exp(0.25)))
    
end
semilogy(N, error,'-o');