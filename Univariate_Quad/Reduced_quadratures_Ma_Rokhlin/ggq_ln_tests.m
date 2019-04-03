clear;clc;close all;

% quad = csvread('ggq_ln_40.csv'); N = 40;
% f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2));
% fint = @(k)(mod(k-1,2)==0).*(2.0/(k+1)) + (mod(k-1,2)==1).*(-1.0/(1 + k/2).^2);
quad = dlmread('ggq_1b3_20.dat')'; N = 20;
f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(x.^(k/2 + 1/3));
fint = @(k) (mod(k-1,2)==0).*(2.0/(k+1)) + (mod(k-1,2)==1).*(3.0/(3*(k/2) + 4.0));
% quad = dlmread('gauss_legendre_15.dat');N=15;
% f = @(x,k) x.^(k-1);
% fint = @(k) 1.0/(k);
x = quad(:,1);w = quad(:,2);


for ii=N:-1:2
   degree = checkDegree(x,w,f,fint)
   [x,w] = reduced_quad(x,w,f); 
   scatter(x,ii*ones(length(x),1),10,w,'filled');
   hold on;
   pause(0.1);
   ylim([0,N+1]);
end