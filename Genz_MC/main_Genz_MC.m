clear;clc;format long;close all;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

nAvg = 10;
d = 3;
for ftype=[1:6]
    numpts = [];
    error = [];
    for i=9:2:20
        eval_pts = rand(i,d);
        avgErr = zeros(nAvg,1);
        for nIter = 1:nAvg
            a = rand(d,1);
            u = rand(d,1);
            [fn, exact] = genz(eval_pts,a,u,ftype);
            mu = mean(fn);
            avgErr(nIter,1) = abs(mu - exact);
        end
        numpts(end+1) = i^d;
        error(end+1) = mean(avgErr);
    end
    figure(ftype)
    loglog(numpts,error,'bo-');
    dlmwrite(strcat('MC_f',num2str(ftype),'.dat'),[numpts',error'],'precision','%2.16f');
end