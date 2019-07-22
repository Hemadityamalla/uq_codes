%Script to demonstrate the error in using piecewise linear approximation
%with linear extrapolation

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');


fn = genz_fns(0:0.001:1,rand(1,1),rand(1,1), 5);
Y = randn(1e6,1);
F = Y(1:10,1);
fpl = griddedInterpolant(sort(F),fn(sort(F)),'linear');
plot(sort(Y),fn(sort(Y)),'b',sort(F),fpl(sort(F)),'ro',sort(Y),fpl(sort(Y)),'g.');
xlabel('\xi'); ylabel('f/f_{pl}'); 
legend('f(\xi)', 'f_{pl}(\xi)', 'f_{pl} over Y');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;