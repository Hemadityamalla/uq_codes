clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

nAvg = 100;
d = 5;
for ftype=[1:5]
    numpts = [];
    error = [];
    for i=5*(10.^[1,2,3,4,5])
        eval_pts = rand(i,d);
        avgErr = zeros(nAvg,1);
        for nIter = 1:nAvg
            a = rand(d,1);
            u = rand(d,1);
            [fn, exact] = genz(eval_pts,a,u,ftype);
            mu = mean(fn);
            avgErr(nIter,1) = abs(mu - exact);
        end
        numpts(end+1) = i;
        error(end+1) = mean(avgErr);
    end
    dlmwrite(strcat('MC_f',num2str(ftype),'.dat'),[numpts',error'],'precision','%2.16f');
end
    
polyfit(log(numpts),log(error),1);