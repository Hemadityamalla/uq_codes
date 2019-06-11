close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 5e2;


%Test function to be integrated
%f = @(x) 1.0*(x > 0.5);
%f = @(x) x.*(x < 0.5) + (x >= 0.5).*exp(-0.5*x);
%f = @(x) exp(-abs(x - 0.5)/2)/sqrt(2*pi);
f = @(x) exp(-x.^2/2/sqrt(2*pi));
error_giq = [];
error_fiq = [];
N = 60;
range = 5:5:N;
D = 20;
numavg = 10;
for degree = range
    degree
    e_giq = 0;
    e_fiq = 0;
    for tt=1:numavg
        ya = rand(Kmax,1);
        [xa,wa,fixedNodes] = cappedaccuracy_quad_v3([D,degree],f,ya);
        [xt,wt,yt] = fixed_implict_quad(degree,Kmax);
        e_giq = e_giq + abs(sum(f(xa).*wa) - mean(f(ya)));
        e_fiq = e_fiq + abs(sum(f(xt).*wt) - mean(f(yt)));
    end
    error_giq = [error_giq;e_giq/numavg];
    error_fiq = [error_fiq;e_fiq/numavg];
end
figure(2);
loglog(range(1:end),error_giq(1:end),'bo-',range(1:end),error_fiq(1:end),'r+-');
xlabel('Number of nodes');ylabel('Error');
legend('adaptive','fixed');
