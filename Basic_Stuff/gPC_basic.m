%clear; clc; format long;

%x = rand(1e4,1);

f = @(x) exp(1 + x);
%gPC approximation
for degree=2
N = degree;
%Generating the quadrature rule
[xi,w] = gaussQuad(3,'Legendre');
%xi = (xi+1.0)/2; w = w/2;
poly_coeffs = gpc_coeffs(N, xi, w, f(xi), 'legendre')
eval_pts = linspace(-1,1,500)';
Y_approx = gpc_polyval(poly_coeffs, eval_pts,'legendre');
plot(eval_pts,f(eval_pts),'r.',eval_pts,Y_approx,'b-');
legend('exact','gPC')
hold on;
grid on
pause(0.5)
end

