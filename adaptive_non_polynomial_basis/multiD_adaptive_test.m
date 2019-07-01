clear;clc; format long;


%dimension 
d = 4;

%points used for piecewise linear interpolation
Kmax = 80;

%Test functions to be integrated
a = rand(d,1);
u = rand(d,1);
ftype = 6;
fn = genz_fns(0:0.001:1, a, u, ftype);

error = [];
degree = 5:2:16;
for N = degree
    %Generating the quadrature rule
    D = 4; S = 5; Y = rand(Kmax,1); 
    [X,W] = multiD_adaptive(fn, N, D, S, Y, d);
    error(end+1) = abs(sum(fn(X).*W) - mean(fn(setprod(Y,d))))

end

semilogy(degree, error);