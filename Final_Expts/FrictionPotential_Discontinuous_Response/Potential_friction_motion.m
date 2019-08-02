%Script solving the stochastic particle moving in a potential field with
%friction

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',0.5,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
%rng(1,'twister');

f = 2; %Larger this value, faster the steady state is reached
N = chebop(@(x, u) diff(u,2) + f*diff(u) + 35*0.5*u.^3 - 15*0.5*u, [0, 10]);



% %This is the Monte-Carlo approach
% xi = rand(5e3,1);
% ic = 0.05 + 0.2*(-1 + 2*xi);
% QoI = [];
% for iter=1:length(ic)
% N.lbc = [ic(iter); 0];
% u = N\0;
% QoI(end+1) = u(end);
% end


%Pseudo-spectral Approach 

[xi, w] = legpts(100); %Degree 19

ic = 0.05 + 0.2*(xi);

QoI = [];
for iter=1:length(ic)
N.lbc = [ic(iter); 0];
u = N\0;
QoI(end+1) = u(end);
end

fhat = gpc_coeffs(3, xi, w', QoI', 'legendre');

%eval_pts = (-0.15:0.01:0.25)'; 
%plot(eval_pts, gpc_polyval(fhat, eval_pts, 'legendre'));
