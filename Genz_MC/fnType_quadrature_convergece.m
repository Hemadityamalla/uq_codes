%Observing the order of convergence of various types of genz functions
%using Gauss-Legendre quadrature

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

addpath('/ufs/hemadity/Documents/chebfun'); %Using chebfun

for fntype = 1:6
   error = [];
   a = rand(1,1);
   u = rand(1,1);
   N = 1:50;
   for nodes = N
      [x,w] = legpts(nodes, [0,1]); %Gauss-legendre quadrature in [0,1] 
      [fn, exact] = genz(x,a,u,fntype);
      error(end+1) = abs(exact - sum(fn.*w'));
   end
   semilogy(N, error,'o-');
   hold on;
end
legend('1','2','3','4','5','6');