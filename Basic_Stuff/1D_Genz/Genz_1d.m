clear;clc;

%Defining the 1D Genz functions

g = @(xi) exp(-0.5*abs(xi - 0.5)); %continuous 
%g = @(xi) 1.0./(0.5^(-2) + (xi - 0.5).^2); %Product Peak 
%g = @(xi) (1 + 0.5*xi).^(-2); %Corner peak 
%g = @(xi) (xi > 0.5)*0 + (xi >= 0.5).*exp(0.5*xi); %Discontinuous
%g = @(xi) exp(-(0.5^2)*(xi - 0.5).^2); %Gaussian peak
%g = @(xi) cos(2*pi*0.5 + 0.5*xi); %Oscillatory

% %Some special G-function
% a = (1 - 2)/2.0;
% g = @(xi) (abs(4*xi - 2) + a)/(1 + a); 




%Gauss-Hermite quadrature points (Q-points)
Q = 50;
[xi,w] = gaussQuad(Q,'Hermite');
g_pts = g(xi); %Integrates polynomials exactly, but not the Genz functions!

ip = (xi(1):0.01:xi(end))'; %Interpolating points

%Number of terms in the polynomial approximation
N = 25;
g_hat = zeros(N,1);
g_approx = 0;
for i=1:N+1
    g_hat(i,1) = sum(w.*hermite(xi,i-1).*g_pts)/(factorial(i-1));
    g_approx = g_approx + g_hat(i,1)*hermite(ip,i-1);
end

plot(xi,g_pts,'b.-');
hold on;
plot(ip, g_approx,'r');
xlim([-5,5]);
ylim([min(g_pts)-0.5,max(g_pts)+0.5]);

