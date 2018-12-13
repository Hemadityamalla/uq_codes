clear;
clc;

f = @(x) cos(2*pi + 5.0*(x(:,1) + x(:,2)));%x(:,1).*exp(x(:,2))./(1 + x(:,3).^2);

d = 2; %dimension of the random vector

polyBasis = 'Hermite';

%Gaussian quadrature points (q^d points)
q = 20;
[xi,w] = gaussQuad(q,polyBasis);
%Change the setprod in such a way that you supply a single array and 'd'-
%DONE
%the unique command in setprod filters out some necessary combinations if
%the initial arrays have the same entries
%DONE
eval_pts = setprod(xi,d);
weights = setprod(w,d);

%Quadrature rule for the mean square error(Q^d points)
Q = 20;
[xi_mse,w_mse] = gaussQuad(Q,polyBasis);
eval_pts_mse = setprod(xi_mse,d);
weights_mse = setprod(w_mse,d);


order = 1;%[0,1,2,3,4,5,6]; %maximum degree of the multivariate polynomial
MSE = [];
for N = order
lexOrdering = monomialDegrees(d,N);
%Pre-computing the normalization factors
gamma = factorial(0:N);
P = size(lexOrdering,1);
fhat = zeros(P,1);
fapprox = 0;
for i_P = 1:P
    fhat(i_P,1) = (sum(hermite(eval_pts,lexOrdering(i_P,:)').*f(eval_pts).*prod(weights,2)))/prod(gamma(lexOrdering(i_P,:)+1));
    fapprox = fapprox + fhat(i_P,1)*hermite(eval_pts_mse, lexOrdering(i_P,:)');
end
%Computing the mean squared error
MSE = [MSE;sqrt(sum((fapprox - f(eval_pts_mse)).^2.*prod(weights_mse,2)))];
end

%plotting, works only for d=2
if d == 2
    [x,y] = meshgrid(xi,xi);
    fn = reshape(f(eval_pts),q,q);
    fn_approx = reshape(fapprox,q,q);
    surf(x,y,fn,'FaceColor',[1,0,0]);
    hold on;
    surf(x,y,fn_approx,'FaceColor',[0,1,0]);
    lim = 2;
    xlim([-lim,lim]);ylim([-lim,lim]);
end