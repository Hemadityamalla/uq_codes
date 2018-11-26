clear;clc;

%Defining the 1D Genz functions

%g = @(xi) exp(-0.5*abs(xi)); %continuous 
%g = @(xi) 1.0./(0.5^(-2) + (xi).^2); %Product Peak 
%g = @(xi) (1 + 5*xi).^(-2); %Corner peak 
g = @(xi) (xi > 0)*0 + (xi <= 0).*exp(0.5*xi); %Discontinuous
%g = @(xi) exp(-(0.5^2)*(xi).^2); %Gaussian peak
g = @(xi) cos(2*pi*0.5 + 5*xi); %Oscillatory

% %Some special G-function
% a = (1 - 2)/2.0;
% g = @(xi) (abs(4*xi - 2) + a)/(1 + a); 




%Gauss-Hermite/Legendre quadrature points (Q-points)
Q = 20;
[xi,w] = gaussQuad(Q,'Legendre');
g_pts = g(xi);

ip = (xi(1):0.01:xi(end))'; %Interpolating points

%Number of terms in the polynomial approximation
N = 19;
g_hat = zeros(N,1);
g_approx = 0;
for i=1:N+1
    polynomial = legendre(xi, i-1);
    g_hat(i,1) = sum(w.*polynomial.*g_pts)/(sum(w.*polynomial.*polynomial));
    g_approx = g_approx + g_hat(i,1)*legendre(ip,i-1);
end

plot(xi,g_pts,'bo-');
hold on;
plot(ip, g_approx);
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(g_pts)-0.5,max(g_pts)+0.5]);

