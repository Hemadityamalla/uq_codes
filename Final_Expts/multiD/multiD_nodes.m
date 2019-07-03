clear;clc; format long;close all;
xpos = 500;ypos = 500; width = 800; height = 800;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);



%points used for piecewise linear interpolation
Kmax = 100;
rng(1,'twister'); %Seeding for reproducibility.

%dimension 
d = 2; %Can't go to 5

for fnType = 1:6
    
    %Test functions to be integrated
    a = rand(d,1);
    u = rand(d,1);
    fn = genz_fns(0:0.001:1, a, u, fnType);
    %Generating the quadrature rule
    N=8;D = 8; S = 5; Y = rand(Kmax,1); 
    [X,W] = multiD_adaptive(fn, N, D, S, Y, d);
    figure(fnType);
    scatter(X(:,1),X(:,2),16,W,'b','filled')%,range(1:end),error_fiq(1:end),'r^-', range(1:end), error_mc(1:end),'g*-');
    xlabel('\xi_1');ylabel('\xi_2');
    %legend('Approx. integrand','Implicit','Monte Carlo');
    grid on;box on;%set(gcf,'Position',[xpos ypos width height]);box on;
end