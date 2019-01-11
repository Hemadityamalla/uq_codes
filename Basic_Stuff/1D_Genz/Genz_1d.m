clear;clc;close all; format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
%Defining the 1D Genz functions

g = @(xi) exp(-0.5*abs(xi)); %continuous 
%g = @(xi) 1.0./(0.5^(-2) + (xi).^2); %Product Peak 
%g = @(xi) (1 + 0.05*xi).^(-2); %Corner peak 
g = @(xi) (xi > 0)*0 + (xi <= 0).*exp(0.5*xi); %Discontinuous
%g = @(xi) exp(-(0.5^2)*(xi).^2); %Gaussian peak
g = @(xi) cos(2*pi*0.5 + 0.5*xi); %Oscillatory

% %Some special G-function
%a = (1 - 2)/2.0;
%g = @(xi) (abs(4*xi - 2) + a)/(1 + a); 


%g = @(xi) abs(xi);

error = [];
error_mean = [];
polyBasis = 'Legendre';
%Gauss-Hermite/Legendre quadrature points (Q-points)
Q = 100;
[xi,w] = gaussQuad(Q,polyBasis);
g_pts = g(xi);


%MSE quadrature points 
Q_mse = 500;
[xi_mse,w_mse] = gaussQuad(Q_mse,polyBasis);


%Number of terms in the polynomial approximation
%order = [1,4,9,14,19,24,29];
order = 1:29;
%dlmwrite('coeff.dat',order);
N_iter = 1;
for N = order
    
    %gamma = factorial(0:N); %for hermite basis
    gamma = 2.0./(2*(0:N) + 1.0); %Legendre basis
    g_hat = zeros(round(N),1); 
    g_approx = 0;
    for i=1:N+1
        polynomial = legendre(xi, i-1);
        %g_hat(i,N_iter) = sum(w.*polynomial.*g_pts)/(sum(w.*polynomial.*polynomial));
        g_hat(i,1) = sum(w.*polynomial.*g_pts)/(gamma(i));
        g_approx = g_approx + g_hat(i,N_iter)*legendre(xi_mse,i-1);
    end
    %dlmwrite('coeff.dat',g_hat','-append');
    
    %error = [error; abs(g_hat(1,1) - gmean)];
   %error = [error; sqrt(sum((exp(-0.5*ip.^2)/sqrt(2*pi)).*abs(g(ip) - g_approx).^2))];
   error_mean = [error_mean;abs((sum(g(xi_mse).*w_mse) - sum(g_approx.*w_mse))/sum(g(xi_mse).*w_mse))];
   error = [error;sqrt(sum(w_mse.*(g(xi_mse) - g_approx).^2))];
hold on;
plot(xi_mse, g_approx,'LineWidth',2);
xlim([min(xi)-0.05,max(xi)+0.05]);
end

plot(xi,g_pts,'bo-');
hold on;
plot(xi_mse, g_approx,'r');
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(g_pts)-0.5,max(g_pts)+0.5]);
figure(2)
loglog(order, error, '-o');
grid on
hold on;
%p = polyfit(log(order),log(error)',1);
%loglog(order,order.^(-0.5*p(1)),'r');
xlabel('N (Number of gPC expansion terms)');ylabel('Mean-square Error');
