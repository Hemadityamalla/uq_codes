clear; clc;
N = 10;
dx = 1.0/N;
x = [0:dx:1]';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
T_left = 1;
c0 = 1.0 ;
% %Define phi- one random variable for SC
% d = 1; 
% xi = randn(d,1);
sigma = 1 ;

%Gauss-Hermite quadrature points (Q points)
Q = 25;
[xi,w] = gaussQuad(Q,'Hermite');


%Deterministic solver

for i=[1:1:Q]
    gmid = xi(i)*((sqrt(2)*sigma)./((0.5)*pi))*(sin(pi*xmid*(0.5)));
    cmid = c0*exp(gmid);
    c_left = cmid(2:N); % note this is because of how Matlab's spdiags works
    c_right = cmid(1:N-1);
    c_diag = -(c_left+c_right);
    A = 1/(dx^2)*spdiags([c_left c_diag c_right],[-1 0 1],N-1,N-1);
    b = zeros(N-1,1) ;
    b(1) = -(cmid(1)/(dx^2))*T_left;
    u = A\b;
    T_half(i,1) = u(N/2);
end

%Evaluating the coefficients

%Hermite polynomials(probabilist)
h0 = @(xi) ones(size(xi));
h1 = @(xi) xi;
h2 = @(xi) xi.^2 - 1.0;
h3 = @(xi) xi.^3 - 3.0*xi;
h4 = @(xi) xi.^4 - 6.0*xi.^2 + 3.0;

T_0 = sum(w.*h0(xi).*T_half)/(factorial(0));
T_1 = sum(w.*h1(xi).*T_half)/(factorial(1));
T_2 = sum(w.*h2(xi).*T_half)/(factorial(2));
T_3 = sum(w.*h3(xi).*T_half)/(factorial(3));