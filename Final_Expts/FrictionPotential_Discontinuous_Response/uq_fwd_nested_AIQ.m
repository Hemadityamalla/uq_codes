clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');


%Loading the MC results
data = dlmread('data_t2.dat');

%Pseudo-spectral Approach
exactmean = mean(data(:,2));
exactstd = std(data(:,2));
err = [];err_std = [];
mu = []; std_sq = []; 
Kmax = 1e6; y = data(:,1);

%Inputs for the nested quadrature rule
lmax = 8;
D = 4; N = 5;

[QR] = AIQ_nested([D,N],lmax,data(:,1)); %Need to store the function evaluations

for ii=1:lmax
   x = QR{ii}.nodes;
   w = QR{ii}.weights;
   QoI = Potential_friction_motion(x);
   mu(end+1) = sum(QoI.*w);
   std_sq(end+1) = sum((QoI - mu(end)).^2.*w);
   err(end+1) = abs(mu(end) - exactmean);
   err_std(end+1) = abs(sqrt(std_sq(end)) - exactstd);
end

