clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',16);

d=3;
nAvg = 10;
for ftype=1:6
    error = [];
    numpts = [];
    for pts= 2:5
        q = (2^pts)-1;
        [xi,w] = clencurt(q);
        xi = 0.5*(xi + 1);
        w = w/2;
        eval_pts = setprod(xi,d);
        weights = prod(setprod(w,d),2);
        avgErr = zeros(nAvg,1);
        for nIter = 1:nAvg
            a = rand(d,1);
            u = rand(d,1);
            [fn, exact] = genz(eval_pts,a,u,ftype);
            soln = dotprod(fn',weights);
            avgErr(nIter,1) = abs(soln - exact);
        end
        numpts(end+1) = length(eval_pts);
        error(end+1) = mean(avgErr);
    end
    figure(ftype)
    loglog(numpts, error,'bo-');
    dlmwrite(strcat('Tensor_f',num2str(ftype),'d',num2str(d),'.dat'),[numpts',error'],'precision','%2.16f');
end
