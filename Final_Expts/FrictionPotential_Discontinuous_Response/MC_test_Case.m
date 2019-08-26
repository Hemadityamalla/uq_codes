clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');


Kmax = 1e5; y = -1 + 2*rand(Kmax,1);
QoI = [];
%Monte Carlo approach for exact value

f = 2; %Larger this value, faster the steady state is reached
N = chebop(@(x, u) diff(u,2) + f*diff(u) + 35*0.5*u.^3 - 15*0.5*u, [0, 2]); %ODE
ic = 0.05 + 0.2*(y);%Initial condition

fevals = zeros(length(y),1);
for iter=1:length(ic)
N.lbc = [ic(iter); 0];
u = N\0;
fevals(iter,1) = u(end);
end
