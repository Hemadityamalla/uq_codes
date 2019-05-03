%clear; clc; format long;

%x = rand(1e4,1);

f = @(x) x;
%gPC approximation
for degree=6
N = degree;
%Generating the quadrature rule
[xi,w] = gaussQuad(10,'Legendre');
%xi = (xi+1.0)/2; w = w/2;
poly_coeffs = gpc_coeffs(N, xi, w, f(xi), 'legendre');
eval_pts = linspace(min(xi),max(xi),50)';
Y_approx = gpc_polyval(poly_coeffs, eval_pts);
% plot(xi,f(xi),'r.',eval_pts,Y_approx,'b-');
% ylim([-0.7,2])
% hold on;
% grid on
% pause(0.5)
end

