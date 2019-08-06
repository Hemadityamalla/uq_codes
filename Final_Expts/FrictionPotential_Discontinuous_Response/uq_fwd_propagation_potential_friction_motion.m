clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');


%Loading the MC results
data = dlmread('MC_fevals.dat');


%Pseudo-spectral Approach
N = 5*[1,2,4,8,16,32];
exactmean = 0.163663;
exactstd = 0.633865691;
err = [];err_std = [];
mu = []; std_sq = []; 
Kmax = 1e6; y = data(1:1e3,1);%-1 + 2*rand(Kmax,1);

D = 10; %Input for the adaptivequad
for i=N
    i
[xi,w] = legpts(i); w = 0.5*w;
%[xi, w] = fixed_implict_quad(i,y);w=w';
QoI = Potential_friction_motion(xi);
%Solving the model as a blackbox
%nodes = y(1:D);L = setdiff(y,nodes);
%QoI = Potential_friction_motion(nodes);

%[xi,w] = fdaqr([i,D], QoI, L, nodes); w = w';
%QoI2 = Potential_friction_motion(setdiff(xi, nodes));

%Postprocessing--- gPC of the output
expansionDegree = 3;

%fhat = gpc_coeffs(expansionDegree, xi, w', [QoI2;QoI], 'legendre');
%QoI = [QoI2;QoI];
mu(end+1) = sum(QoI.*(w)');
%gamma = 2.0./(2*(0:expansionDegree)' + 1.0);
std_sq(end+1) = sum((QoI - mu(end)).^2.*(w)');
err(end+1) = abs(exactmean - mu(end));
err_std(end+1) = abs(exactstd - sqrt(std_sq(end)));
end
loglog(N,err,'r-o',N,err_std,'r-*');
