%Script to test the convergence of f_pl as it's accuracy is increased

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');

fn = genz_fns(0:0.001:1,rand(1,1),rand(1,1), 6);
Y = rand(1e6,1);
N = 2.^[1:5];
error = [];
for D = N
    F = 0:(1./D):1;%Y(1:D,1);
    evals = min(F):0.001:max(F);
    fpl = griddedInterpolant(sort(F),fn(sort(F)'),'linear');
    plot(evals',fn(evals'),'b',evals,fpl(evals),'ro')%,sort(Y),fpl(sort(Y)),'g.')
    error(end+1) = norm(fn(evals') - fpl(evals)');
end
loglog(N, error,'bo-',N,1./N,'r-');
%xlabel('\xi'); ylabel('f/f_{pl}'); 
%legend('f(\xi)', 'f_{pl}(\xi)', 'f_{pl} over Y');
%grid on;set(gcf,'Position',[xpos ypos width height]); box on;

polyfit(log(1./N),log(error),1)