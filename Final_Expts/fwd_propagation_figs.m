clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

xpos = 500;ypos = 500; width = 1000; height = 800;

addpath('/ufs/hemadity/Documents/chebfun'); %Using chebfun


%Output

x = 0:0.1:20;

k= 2; l=5;

output = (k/l)*(x/l).^(k-1).*exp(-(x/l).^k) ;%+ 0.01*rand(1,length(x));
figure(1)
area(x,output);

%Inputs

x = -10:0.01:10;

input1 = (1/sqrt(2*pi))*exp(-0.5*x.^2);
figure(2)
area(x,input1);

input2 = (1/(sqrt(2*3*pi)))*exp(-(x.^2)/(2*3));
figure(3)
area(x,input2)
ylim([0, 0.4])