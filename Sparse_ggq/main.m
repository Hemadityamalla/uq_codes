clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',16);

d=2;
nAvg = 100;
quadrature = 'ggq';
for ftype=[1:6]
    error = [];
    numpts = [];
    for pts=[1,2,3,4,5]
        growth = @(x) 2^(x-1);
        [xi,w] = smolyakSparseGrid(d,pts,growth, quadrature);
        if strcmp(quadrature,'ggq')
            eval_pts = xi'; weights = w';
        else
            eval_pts = 0.5*(1.0 + xi'); weights = w'/2^d;
        end
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
    dlmwrite(strcat('Smolyak_',quadrature,'_f',num2str(ftype),'d',num2str(d),'.dat'),[numpts',error'],'precision','%2.16f');
end
