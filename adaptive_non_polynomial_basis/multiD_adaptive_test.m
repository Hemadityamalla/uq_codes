clear;clc; format long;close all;


%dimension 
d = 2;

%points used for piecewise linear interpolation
Kmax = 80;

%Test functions to be integrated
a = rand(d,1);
u = rand(d,1);
ftype = 5;
fn = genz_fns(0:0.001:1, a, u, ftype);

error = [];
degree = 5:2:25;
numavg =  10;
for N=degree
    err = 0;
    N
    for i = 1:numavg
        %Generating the quadrature rule
        D = 4; S = 5; Y = rand(Kmax,1); 
        [X,W] = multiD_adaptive(fn, N, D, S, Y, d);
        err = err + abs(sum(fn(X).*W) - mean(fn(setprod(Y,d))));

    end
    error(end+1) = err/numavg;
end

semilogy(degree, error,'bo-');