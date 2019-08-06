clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
%addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;

data = dlmread('MC_fevals.dat');

GL = dlmread('GL_error.dat');
FIQ = dlmread('FIQ_error.dat');
AIQ = dlmread('AIQ_error.dat');
MC = dlmread('MC_error.dat');

exactmean = 0.163663;
exactstd = 0.633865691;

N = 5*[1,2,4,8,16,32];
err = []; err_std = [];
for i=N
    err(end+1) = abs(exactmean - mean(data(1:i,2)));
    err_std(end+1) = abs(exactstd - std(data(1:i,2)));
end

%Mean error
figure(1)
loglog(GL(:,1),GL(:,2),'ro-',N, err, 'bo-', AIQ(:,1),AIQ(:,2),'go-');
hold on;grid on;set(gcf,'Position',[xpos ypos width height]); box on;
xlabel('N'); ylabel('Abs. error');legend('Gauss-Legendre','Monte Carlo','Approximate-integrand');

%sigma error
figure(2)
loglog(GL(:,1),GL(:,3),'r*-',N,err_std,'b*-',AIQ(:,1),AIQ(:,3),'g*-');
hold on;grid on;set(gcf,'Position',[xpos ypos width height]); box on;
xlabel('N'); ylabel('Abs. error');legend('Gauss-Legendre','Monte Carlo','Approximate-integrand');