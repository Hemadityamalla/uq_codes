clear;clc; format long;close all;
xpos = 500;ypos = 500; width = 1000; height = 800;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);



%points used for piecewise linear interpolation
Kmax = 80;
rng(1,'twister'); %Seeding for reproducibility.

%dimension 
d = 5; %Can't go to 5

for fnType = 1:6
    
    %Test functions to be integrated
    a = rand(d,1);
    u = rand(d,1);
    testpts = rand(10,d);
    fn = genz_fns(testpts, a, u, fnType);
    error = [];
    range = 5:2:30;
    numavg =  50;
    for N=range
        err = 0;
        N
        for i = 1:numavg
            %Generating the quadrature rule
            D = 4; S = 5; Y = rand(Kmax,1); 
            [X,W] = multiD_adaptive(fn, N, D, S, Y, d);
            err = err + abs(sum(fn(X).*W) - mean(fn(setprod(Y,d))));
            %err = err + abs(sum(fn(X).*W) - mean(fn(rand(1e4,d)))); %This is incorrect, hence does not give convergence

        end
        error(end+1) = err/numavg;
    end
    %Loading data
    data_MC = dlmread(strcat('MC_f',num2str(fnType),'.dat'));
    data_tensor = dlmread(strcat('Tensor_f',num2str(fnType),'d',num2str(d),'.dat'));
    figure(fnType);
    loglog(range(1:end).^d,error(1:end),'bo-',data_MC(1:end,1),data_MC(1:end,2),'r^-', data_tensor(1:end,1), data_tensor(1:end,2),'g*-');
    xlabel('Number of nodes');ylabel('Error');
    legend('Approx. integrand','Monte Carlo','Tensor grid');
    grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end