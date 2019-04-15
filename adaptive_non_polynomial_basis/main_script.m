clear;clc; format long;


%points used for piecewise linear interpolation
Kmax = 5000;



f = @(x) exp(-x.^2/2)/sqrt(2*pi);%Test function to be integrated

exact = 0.5*erf(1.0/sqrt(2));
error = [];
for degree = 2:20
    [x,w] = adaptive_quad(degree,f,Kmax);
    error = [error;abs(dotprod(f(x'),w)-exact)];
end

plot(2:20,error);
