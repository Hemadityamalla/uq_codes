clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);


testFn = 6;


Kmax = 2.5e2;
error = [];
degree = [2:2:30];
for N = degree
    N
    %log(x)
    %f_deets.fn = @(x,k) (k == 1).*(x.^(k-1)) + (k==2).*(log(x)) +(k > 2).*(x.^(k-1).*log(x)); f_deets.coeffs = 1:N;
    %sin(x)
    f_deets.fn = @(x,k) (k == 1).*(x.^(k-1)) + (k==2).*(sin(x)) +(k > 2).*(x.^(k-1).*sin(x)); f_deets.coeffs = 1:N;
    %x^k
    %f_deets.fn = @(x,k) x.^(k-1); f_deets.coeffs = 1:N; 
    %haar
    %f_deets.fn = @(x,k) (k==1).*x.^(k-1) + (k > 1).*(haar(k-1,x));f_deets.coeffs = 1:N;
    %f(polynomials)
    %f_deets.fn = @(x,k) (k == 1).*(x.^(k-1)) + (k > 1).*(log(x.^(k-1))); f_deets.coeffs = 1:N;
    mean_error = 0.0;
    Navg = 50;
    for iter = 1:Navg
        %Generalized Implicit Quadrature rule
        [x,w,y] = general_fixed_implicit_quad(f_deets, N, Kmax);
        %Integrand
        a = rand(1); u = rand(1);
        [fn_quad, ~] = genz(x,a,u, testFn);
        [fn_samples, ~] = genz(y,a,u, testFn);
        exact = mean(fn_samples);
        quadrature = sum(fn_quad.*w);
        mean_error = mean_error + abs(quadrature - exact);
    end
    mean_error = mean_error/Navg;
    error(end+1) = mean_error;
end
errpoly = dlmread('poly_genz1.dat');
loglog(degree, error,'bo-',errpoly(:,1),errpoly(:,2),'r*-');
