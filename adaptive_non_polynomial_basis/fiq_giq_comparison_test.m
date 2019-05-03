%Script to check if supplying information about the function to the QR
%generator helps increasing the order of convergence

clear;clc; format long;

%Sample points
Kmax = 1e3;

interp_pts = linspace(0,1,25);
fn = @(x) exp(-x.^2/2/sqrt(2*pi));%(x.*(x < 0.5) + (x >= 0.5).*exp(-5*x)); 
pp = griddedInterpolant(interp_pts, fn(interp_pts),'linear');
f_deets.fn = @(x,k) (k==2)*fn(x) + (k==1).*1 + (k>2).*x.^(k-2);
f_deets.coeffs = [];
fnmarker.fnidx = 2; fnmarker.pts = interp_pts;

error_giq = [];
error_fiq = [];
D = 75;
range = 5:5:D;
numavg = 10;
for degree = range
    degree
    e_giq = 0;
    e_fiq = 0;
    for tt=1:numavg
    [xa,wa,ya] = approxgeneral_fixed_implicit_quad(f_deets,fnmarker,degree,Kmax);
    [xt,wt,yt] = fixed_implict_quad(degree,Kmax);
    e_giq = e_giq + abs(dotprod(f_deets.fn(xa',2),wa)-mean(f_deets.fn(ya,2)));
    e_fiq = e_fiq + abs(mean(f_deets.fn(yt,2)) - dotprod(f_deets.fn(xt',2),wt));
    end
    error_giq = [error_giq;e_giq/numavg];
    error_fiq = [error_fiq;e_fiq/numavg];
end
%figure(2);
%semilogy(range,error_giq,'b',range,error_fiq,'r');
figure(1);
xx = 0:.001:1;
plot(xx, fn(xx));
hold on
plot(interp_pts, pp(interp_pts), 'o-');

figure(2);
loglog(range(1:end),error_giq(1:end),'bo-',range(1:end),error_fiq(1:end),'r+-');
xlabel('Number of nodes');ylabel('Error');
legend('adaptive','fixed');

%Expt. 1
%Supplying the function into the Vandermonde matrix obviously gives a QR
%that exactly integrates the function. So the error in such a case is
%obviously zero for all degrees. Two cases have been tested: 1) A C_inf
%Gaussian function, 2) A discontinuous step function. 
%Obs. 1
%Both the cases yield the same result as described in the second sentence.

%Expt. 2
%In most real cases, the 'exact function' is not available. So an
%approximation of the function is supplied as a basis to the QR generator.
%The type of approximation used and its effect on the convergence of the QR
% is a matter of further investigation. As a first step, a piecewise linear
% approximation in equispaced points is taken. Then it has been planned to
% take spline, pchip, etc. 
%Obs. 2 
%There appears to be two effects at play here- type of function being
%approximated, accuracy of the approximation. These two effects may or may
%not be `correlated'. Further investigation needed.

%Expt. 3
%Three types of functions are used: Gaussian (mean at 0.5), C_0 fn with peak at 0.5, disc. function at 0.5.
%With increasing order of accuracy (h=0.25,0.1,0.05) of piecewise linear approximation, their
%errors are investigated (based on the observations of Expt. 2).
%Obs. 3
%C_inf function- FIQ is waaaaaay better
%C0 fn.- AIQ is better by a small order(data needs to be averaged).
%Disc. fn.- AIQ performs similar to that of FIQ (data needs to be averaged).
% What happens if I replace a different row or what happens if I replace
% more than one row and with what kinda functions?


