clear;clc;
%Defining the 1D Genz functions

g = @(xi) exp(abs(xi)); %continuous 
%g = @(xi) 1.0./(0.5^(-2) + (xi).^2); %Product Peak 
%g = @(xi) (1 + 5*xi).^(-2); %Corner peak 
%g = @(xi) (xi > 0)*0 + (xi <= 0).*exp(0.5*xi); %Discontinuous
%g = @(xi) exp(-(0.5^2)*(xi).^2); %Gaussian peak
%g = @(xi) cos(2*pi*0.5 + 5*xi); %Oscillatory

% %Some special G-function
% a = (1 - 2)/2.0;
% g = @(xi) (abs(4*xi - 2) + a)/(1 + a); 


%g = @(xi) xi.^4;


%Gauss-Hermite/Legendre quadrature points (Q-points)
Q = 50;
[xi,w] = gaussQuad(Q,'Legendre');
g_pts = g(xi);

ip = (xi(1):0.01:xi(end))'; %Interpolating points

%Number of terms in the polynomial approximation
N = 2;
error = [];
for N = [1:2:20]
g_hat = zeros(round(N),1);
g_approx = 0;
for i=1:N+1
    polynomial = legendre(xi, i-1);
    g_hat(i,1) = sum(w.*polynomial.*g_pts)/(sum(w.*polynomial.*polynomial));
    g_approx = g_approx + g_hat(i,1)*legendre(ip,i-1);
end
%error = [error; sqrt(sum((exp(-0.5*ip.^2)/sqrt(2*pi)).*abs(g(ip) - g_approx).^2))];
end

plot(xi,g_pts,'bo-','LineWidth',2);
hold on;
plot(ip, g_approx,'r','LineWidth',2);
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(g_pts)-0.5,max(g_pts)+0.5]);
%figure(2)
%semilogy([1:2:20], error, '-o');