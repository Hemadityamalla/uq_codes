%Observing the order of convergence of various types of genz functions
%using Gauss-Legendre quadrature

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

xpos = 500;ypos = 500; width = 1000; height = 800;

addpath('/ufs/hemadity/Documents/chebfun'); %Using chebfun

for fntype = 1:6
   error = [];
   a = rand(1,1);
   u = rand(1,1);
   N = 1:50;
   for nodes = N
      [x,w] = legpts(nodes, [0,1]); %Gauss-legendre quadrature in [0,1] 
      [fn, exact] = genz(x,a,u,fntype);
      error(end+1) = abs(exact - sum(fn.*w'));
   end
   loglog(N, error,'o-');
   hold on;
end
legend('f_1','f_2','f_3','f_4','f_5','f_6');ylabel('Absolute Error');xlabel('Nodes')
grid on;set(gcf,'Position',[xpos ypos width height]); box on;