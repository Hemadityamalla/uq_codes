clear;clc;close all;

addpath('/ufs/hemadity/Documents/chebfun');

N = 4; d = 2;
polyspan = monomialDegrees(d,N); %Span size= (N+d)!/(N!*d!)

[x,w] = chebpts(N);

%Multivariate Vandermonde matrix



