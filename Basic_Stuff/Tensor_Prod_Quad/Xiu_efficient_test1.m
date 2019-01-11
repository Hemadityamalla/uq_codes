clear;clc;format long;

f = @(x) exp(-0.5*(abs(x(:,1)) + abs(x(:,2))));%x(:,1).*exp(x(:,2))./(1 + x(:,3).^2);
f = @(x) cos(2*pi*0.5 + 0.5*(x(:,1) + x(:,2)));
f = @(x) (x(:,1) > 0)*0 .* (x(:,2) > 0)*0 + (x(:,1) <= 0).*exp(0.5*(x(:,1) + x(:,2))).*(x(:,2) <= 0);
d = 2; %dimension of the random vector
% Input for type of orthogonal polynomial basis. Alternatively, use 'Hermite'
polyBasis = 'Legendre';

q = 10;
[xi,w] = gaussQuad(q,polyBasis);
eval_pts = setprod(xi,d); %(q^d points)
weights = setprod(w,d);

%Quadrature rule for the mean square error(Q^d points)
Q = 50;
[xi_mse,w_mse] = gaussQuad(Q,polyBasis);
eval_pts_mse = setprod(xi_mse,d);
weights_mse = setprod(w_mse,d);

order = [0,1,2,3,4,5]; %maximum degree of the multivariate polynomial
MSE = []; %Empty array to store the mean-squared error
for N = order
     
    lexOrdering = monomialDegrees(d,N);
    %Pre-computing the normalization factors
    if strcmp(polyBasis,'Hermite')
        gamma = factorial(0:N); %Hermite
    else
        gamma = 2.0./(2*(0:N) + 1.0); %Legendre
    end
    P = size(lexOrdering,1);
    fhat = zeros(P,1);
    fapprox = 0;
    for i_P = 1:P
        fhat(i_P,1) = (sum(legendre(eval_pts,lexOrdering(i_P,:)').*f(eval_pts).*prod(weights,2)))/prod(gamma(lexOrdering(i_P,:)+1));
        fapprox = fapprox + fhat(i_P,1)*legendre(eval_pts_mse, lexOrdering(i_P,:)');
    end
    %Computing the mean squared error
    MSE = [MSE;sqrt(sum((fapprox - f(eval_pts_mse)).^2.*prod(weights_mse,2)))];
end
% This part is just for visualization, works only for d=2
if d == 2
    [x,y] = meshgrid(xi,xi);
    fn = reshape(f(eval_pts),q,q);
    fn_approx = reshape(fapprox,Q,Q);
    surf(x,y,fn,'FaceColor',[1,0,0]); %Exact function is in red
    hold on;
    [x_a,y_a] = meshgrid(xi_mse,xi_mse);
    surf(x_a,y_a,fn_approx,'FaceColor',[0,1,0]); %Approx. function is in green
    lim = 2;
    xlim([-lim,lim]);ylim([-lim,lim]);zlim([-lim,lim])
end
