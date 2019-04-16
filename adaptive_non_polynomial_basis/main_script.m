clear;clc; format long;


%points used for piecewise linear interpolation
Kmax = 1000;

f = @(x) 1.0*(x > 0.5);%Test function to be integrated

exact = 0.5;
error_adaptive = [];
error_trad = [];
for degree = 5:5:100
    [xa,wa] = adaptive_quad(degree,f,Kmax);
    [xt,wt] = fixed_implict_quad(degree,Kmax);
    error_adaptive = [error_adaptive;abs(dotprod(f(xa'),wa)-exact)];
    error_trad = [error_trad;abs(exact - dotprod(f(xt'),wt))];
end
figure(2);
semilogy(5:5:60,error_adaptive,'b',5:5:60,error_trad,'r');
