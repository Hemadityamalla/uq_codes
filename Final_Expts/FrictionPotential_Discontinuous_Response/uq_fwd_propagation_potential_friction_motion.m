clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
%addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');


%Loading the MC results
data = dlmread('MC_fevals.dat');


%Pseudo-spectral Approach
N = 5*[1,2,4,8,16,32,64,128,256];
exactmean = 0.16366;
err = [];
mu = []; std_sq = []; 
Kmax = 1e6; y = data(:,1);%-1 + 2*rand(Kmax,1);

D = 10; %Input for the adaptivequad
for i=N
[xi, w] = fixed_implict_quad(i,y);w=w';
%Solving the model as a blackbox
% nodes = y(1:D);L = setdiff(y,nodes);
% QoI = Potential_friction_motion(nodes);
% 
% [xi,w] = fdaqr([i,D], QoI, L, nodes); w = w';
% QoI2 = Potential_friction_motion(setdiff(xi, nodes));
%[xi,w] = legpts(i);
QoI = Potential_friction_motion(xi);

%Postprocessing--- gPC of the output
expansionDegree = 3;
fhat = gpc_coeffs(expansionDegree, xi, w', QoI, 'legendre');
mu(end+1) = fhat(1);
gamma = 2.0./(2*(0:expansionDegree)' + 1.0);
std_sq(end+1) = sum(gamma.*fhat.^2) - gamma(1)*fhat(1)^2;
err(end+1) = abs(exactmean - fhat(1));
end
%loglog(N,err,'-o');
