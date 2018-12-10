clear;
clc;

f = @(x) x(:,1);%x(:,1).*exp(x(:,2))./(1 + x(:,3).^2);

d = 2; %dimension of the random vector

%Gauss-Hermite quadrature points (Q = q^d points)
q = 8;
[xi,w] = gaussQuad(q,'Legendre');
%Change the setprod in such a way that you supply a single array and 'd'-
%the unique command in setprod filters out some necessary combinations if
%the initial arrays have the same entries
eval_pts = setprod(xi,d);
weights = setprod(w,d);


N = 3; %maximum degree of the multivariate polynomial
lexOrdering = monomialDegrees(d,N);
%Pre-computing the normalization factors
gamma = factorial(0:N);


P = length(lexOrdering);
fhat = zeros(P,1);
fapprox = 0;
for i_P = 1:P
    fhat(i_P,1) = (sum(legendre(eval_pts,lexOrdering(i_P,:)').*f(eval_pts).*prod(weights,2)))/prod(gamma(lexOrdering(i_P,:)+1));
    fapprox = fapprox + fhat(i_P,1)*legendre(eval_pts, lexOrdering(i_P,:)');
end

%plotting(works only for d=2)
[x,y] = meshgrid(xi,xi);
fn = reshape(f(eval_pts),q,q);
fn_approx = reshape(fapprox,q,q);
surf(x,y,fn,'FaceColor',[1,0,0]);
hold on;
surf(x,y,fn_approx,'FaceColor',[0,1,0]);
lim = 2;
xlim([-lim,lim]);ylim([-lim,lim]);