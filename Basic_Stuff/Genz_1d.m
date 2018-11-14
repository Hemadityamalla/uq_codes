clear;clc;

%Defining the 1D Genz functions

g_continuous = @(xi) exp(-0.5*abs(xi - 0.5));

%Gauss-Hermite quadrature points (Q points)
Q = 100;
[xi,w] = gaussQuad(Q,'Hermite');
g_pts = g_continuous(xi);

ip = (xi(1):0.01:xi(end))'; %Interpolating points

%Number of terms in the polynomial approximation
N = 15;
g_hat = zeros(N,1);
g_approx = 0;
for i=1:N
    g_hat(i,1) = sum(w.*hermite(xi,i).*g_pts)/(factorial(i));
    g_approx = g_approx + g_hat(i,1)*hermite(ip,i);
end

plot(xi,g_pts,'.-');
hold on;
plot(ip, g_approx,'r');

