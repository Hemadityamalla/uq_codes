clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
%addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;

data = dlmread('MC_fevals.dat');

GL = dlmread('GL_error.dat');
FIQ = dlmread('FIQ_error.dat');
AIQ = dlmread('AIQ_error.dat');
MC = dlmread('MC_error.dat');

%Mean error
figure(1)
loglog(GL(:,1),GL(:,2),'ro-',FIQ(:,1),FIQ(:,2),'bo-',AIQ(:,1),AIQ(:,2),'yo-');

%sigma error
figure(2)
loglog(GL(:,1),GL(:,3),'ro-',FIQ(:,1),FIQ(:,3),'bo-',AIQ(:,1),AIQ(:,3),'yo-');