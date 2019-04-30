clear;clc; format long;


%points used for piecewise linear interpolation
Kmax = 5e3;

%f = @(x) 1.0*(x > 0.5);%Test function to be integrated
%f = @(x) x.*(x < 0.5) + (x >= 0.5).*exp(-0.5*x);
f = @(x) exp(-abs(x - 0.5)/2)/sqrt(2*pi);

error_adaptive = [];
error_trad = [];
D = 150;
range = 5:10:D;
for degree = range
    [xa,wa,ya] = cappedaccuracy_quad([degree, 10],f,Kmax);
    [xt,wt,yt] = fixed_implict_quad(degree,Kmax);
    error_adaptive = [error_adaptive;abs(dotprod(f(xa'),wa)-mean(f(ya)))];
    error_trad = [error_trad;abs(mean(f(yt)) - dotprod(f(xt'),wt))];
end
figure(2);
semilogy(range,error_adaptive,'b',range,error_trad,'r');
