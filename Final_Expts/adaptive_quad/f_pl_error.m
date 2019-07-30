%Script to demonstrate the error in using piecewise linear approximation
%with linear extrapolation

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');


fn = genz_fns(0:0.001:1,rand(1,1),rand(1,1), 5);
Y = randn(1e6,1);
D = 10;
F = Y(1:D,1);
fpl = griddedInterpolant(sort(F)',fn(sort(F)),'linear');
%Constructing lagrange interpolant (polynomial of degree D)
f_lagrange = chebfun.interp1(F, fn(F));
plot(sort(Y),fn(sort(Y)),'b', sort(Y), f_lagrange(sort(Y)),'k-');
xlabel('\xi'); ylabel('f/f_{lagrange}'); ylim([min(fn(sort(Y))),1]);
legend('f(\xi)', 'f_{lagrange}(\xi)');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;