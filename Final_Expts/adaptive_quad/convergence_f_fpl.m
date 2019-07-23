%Script comparing the convergence rates of integration of a function and
%its piecewise linear interpolant


clear;clc;close all;format long;

addpath('/ufs/hemadity/Documents/chebfun'); %Loading chebfun

%Samples
Y = randn(1e6,1);
F = Y(1:40,1);


fn = genz_fns(0:0.001:1,rand(1,1),rand(1,1), 5); %Test function

fpl = griddedInterpolant(sort(F),fn(sort(F)),'linear'); %Function interpolant


%Computing the exact integral using chebfun
exact = sum(chebfun(fn,[0,1]));

N = 1:50;
e_fn = []; e_pl = [];
for degree=N
   [x,w] = legpts(degree,[0,1]);
   e_fn(end+1) = abs(exact - sum(fn(x).*w'));
   e_pl(end+1) = abs(exact - sum(fpl(x).*w'));
end

loglog(N, e_fn,'bo-',N,e_pl,'ro-');
xlabel('Degree');ylabel('Abs. error');legend('f(x)','f_{pl}(x)');
