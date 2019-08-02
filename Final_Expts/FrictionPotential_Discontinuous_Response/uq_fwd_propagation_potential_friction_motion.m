clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');

%Pseudo-spectral Approach
N = 10:5:30;
exactmean = 0.16366;
err = [];
Kmax = 500; y = -1 + 2*rand(Kmax,1);
for i=N
[xi, w] = fixed_implict_quad(i,y);w=w';
%Solving the model as a blackbox
QoI = Potential_friction_motion(xi);

%Postprocessing--- gPC of the output
fhat = gpc_coeffs(1, xi, 2*w', QoI, 'legendre');
err(end+1) = abs(exactmean - fhat(1));
end
loglog(N,err,'-o');
