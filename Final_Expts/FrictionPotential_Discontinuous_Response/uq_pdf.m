clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/home/hemaditya/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');

data = dlmread('MC_fevals.dat');

exact_soln = @(x) (x > -0.25)*(sqrt(15/35)) - (x < -0.25)*(sqrt(15/35));

%Pseudo-spectral Approach
N = 20; D = 10; %Input for the adaptivequad
y = data(:,1);
nodes = y(1:D);L = setdiff(y,nodes);
QoI = Potential_friction_motion(nodes);

[xi,w] = fdaqr([N,D], QoI, L, nodes); w = w';
QoI2 = Potential_friction_motion(setdiff(xi, nodes));

QoI = [QoI2;QoI];

fhat = gpc_coeffs(3, xi, w', QoI, 'legendre');
samples = -1 + 2*rand(1e7,1);
Xevals = gpc_polyval(fhat, samples, 'legendre');
figure(1)
x = (-1:0.001:1)'; 
plot(x,gpc_polyval(fhat, x, 'legendre'),'b',x,exact_soln(x),'r');
xlim([-1, 1]); xlabel('\xi'); ylabel('X(t\rightarrow \infty, \xi)');
hold on;grid on;set(gcf,'Position',[xpos ypos width height]); box on;legend('gPC','Exact');
figure(2)
histogram(Xevals,'Normalization','probability');
xlim([-1, 1]); xlabel('\xi'); ylabel('PDF(X)');
hold on;grid on;set(gcf,'Position',[xpos ypos width height]); box on;
histogram(exact_soln(samples),'Normalization','probability');
legend('gPC','Exact');



