close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 2.5e2;

for testFn = 1:6
    error_giq = [];
    error_fiq = [];
    N = 25;
    range = 5:5:N;
    D = 5;
    numavg = 10;
    for degree = range
        degree
        e_giq = 0;
        e_fiq = 0;
        for tt=1:numavg
            ya = rand(Kmax,1);
            %Generating the function
            a = rand(1); u = rand(1); 
            f = genz_fns(0:0.01:1, a, u, testFn);
            [xa,wa,fixedNodes] = cappedaccuracy_quad_v3([D,degree],f,ya);
            [xt,wt] = fixed_implict_quad(degree,ya);
            e_giq = e_giq + abs(sum(f(xa).*wa) - mean(f(ya)));
            e_fiq = e_fiq + abs(sum(f(xt).*wt) - mean(f(ya)));
        end
        error_giq = [error_giq;e_giq/numavg];
        error_fiq = [error_fiq;e_fiq/numavg];
    end
    
    figure(testFn);
    loglog(range(1:end),error_giq(1:end),'bo-',range(1:end),error_fiq(1:end),'r+-');
    xlabel('Number of nodes');ylabel('Error');
    legend('adaptive','fixed');
end