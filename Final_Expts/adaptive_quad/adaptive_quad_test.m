close all;clear;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

%points used for piecewise linear interpolation
Kmax = 2.5e2;
rng(1,'twister'); %Seeding for reproducibility.


xpos = 500;ypos = 500; width = 800; height = 800;

for testFn = 5
    error_giq = [];
    ngiq = [];
    error_fiq = [];
    error_mc = [];
    N = 35;
    
    D = 12;
    range = 5:5:N;
    numavg = 15;
    for degree = range
        degree
        e_giq = 0;
        e_fiq = 0;
        e_mc = 0;
        for tt=1:numavg
            ya = rand(Kmax,1);
            %Generating the function
            a = rand(1); u = rand(1); 
            f = genz_fns((0:0.01:1), a, u, testFn);
            [xa,wa,fixedNodes] = cappedaccuracy_quad_v3([D,degree],f,ya);
            [xt,wt] = fixed_implict_quad(length(xa)-1,ya); %degree-1 is used because the approx quadrature rule take degree-1 polynomials as well
            e_giq = e_giq + abs(sum(f(xa).*wa) - mean(f(ya)));
            e_fiq = e_fiq + abs(sum(f(xt).*wt) - mean(f(ya)));
            %Monte-Carlo estimate
            e_mc = e_mc + abs(mean(f(rand(length(xa),1)) - mean(f(ya))));
        end
        ngiq(end+1) = length(xa);
        error_giq = [error_giq;e_giq/numavg];
        error_fiq = [error_fiq;e_fiq/numavg];
        error_mc = [error_mc; e_mc/numavg];
    end
    
    figure(testFn);
    loglog(ngiq,error_giq(1:end),'bo-',range(1:end),error_fiq(1:end),'r^-', range(1:end), error_mc(1:end),'g*-');
    xlabel('Number of nodes');ylabel('Error');
    legend('Approx. integrand','Implicit','Monte Carlo');
    grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end