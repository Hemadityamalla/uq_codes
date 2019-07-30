%Script solving the stochastic, steady state burgers' equation

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',0.5,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');

d = [-1, 1];
x = chebfun('x',d);
N = chebop(d);

N.lbc = 1; N.rbc = -1;

v = 0.1 + exp(log(0.05) + randn(1e3,1)*((log(10))/(2.85)));
QoI = []; wmean = 0; wtotal = zeros(length(-1:0.01:1),length(v));
for iter = 1:length(v)
    N.op = @(w) w.*diff(w) - v(iter)*diff(w,2);
    rhs = 0.0;
    w = N\rhs;
    wtotal(:,iter) = w(-1:0.01:1); 
    QoI(end+1) = w(0.5);
    %plot(w);
    %hold on;
end
wmean = mean(wtotal,2);
plot(-1:0.01:1,wmean,'k');
hold on;
SolVar = var(wtotal')';
curve1 = wmean + 2*sqrt(SolVar);
curve2 = wmean - 2*sqrt(SolVar);
plot(-1:0.01:1,curve1,'b');
hold on;
plot(-1:0.01:1, curve2,'b');
hold on;
%patch([-1:0.01:1, fliplr(curve1)'], [-1:0.01:1, fliplr(curve2)'], 'b');
x2 = [-1:0.01:1, fliplr(-1:0.01:1)];
inBetween = [curve1', fliplr(curve2')];
fill(x2', inBetween','b','FaceAlpha',0.3);
