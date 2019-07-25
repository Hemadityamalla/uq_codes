%Script solving the stochastic, steady state burgers' equation

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');

d = [-1, 1];
x = chebfun('x',d);
N = chebop(d);

N.lbc = 1; N.rbc = -1;

v = 0.05 + exp(log(0.05) + randn(5,1)*((log(10))/(2.85)));
QoI = []; wmean = 0;
for nu = v'
    N.op = @(w) w.*diff(w) - nu*diff(w,2);
    rhs = 0.0;
    w = N\rhs;
    wmean = wmean + w(-1:0.1:1);
    QoI(end+1) = w(0.5);
    %plot(w);
    %hold on;
end

plot(wmean/5);
