close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 5e2;

f = @(x) 1.0*(x > 0.5);%Test function to be integrated
%f = @(x) x.*(x < 0.5) + (x >= 0.5).*exp(-0.5*x);
%f = @(x) exp(-abs(x - 0.5)/2)/sqrt(2*pi);
%f = @(x) exp(-x.^2/2/sqrt(2*pi));

lmax = 15;
D = 5; N = D;
numavg = 25;
efiqavg = zeros(lmax,1)'; egiqavg = zeros(lmax,1)';
for iavg = 1:numavg
    error_giq = [];
    error_fiq = [];
    %Generating nested QR sequence 
    [QR,ya] = cappedaccuracy_quad_nested([D,N], lmax,f, Kmax);
    for ii=1:lmax
        %Generating the FIQR (Gives a QR of degree similar to the ii th QR of the nested sequence)
        [xf,wf,yf] = fixed_implict_quad((N-2)+ii, Kmax);
        error_fiq(end+1) = abs(sum(f(xf).*wf) - mean(f(yf)));
        xa = QR{ii}.nodes; wa = QR{ii}.weights;
        error_giq(end+1) = abs(sum(f(xa).*wa) - mean(f(ya)));
    end
    efiqavg = efiqavg + error_fiq; egiqavg = egiqavg + error_giq;
end
figure(2);
loglog(1:lmax,egiqavg/numavg,'bo-',1:lmax,efiqavg/numavg,'r+-');
xlabel('Number of nodes/levels');ylabel('Error');
legend('adaptive','fixed');
