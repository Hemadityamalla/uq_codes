clear;clc;close all; format long;

%Defining the 1D Genz functions

g = @(xi) exp(-0.5*abs(xi)); %continuous 
%g = @(xi) 1.0./(0.01^(-2) + (xi).^2); %Product Peak 
%g = @(xi) (1 + 0.05*xi).^(-2); %Corner peak 
%g = @(xi) (xi > 0)*0 + (xi <= 0).*exp(0.5*xi); %Discontinuous
%g = @(xi) exp(-(0.5^2)*(xi).^2); %Gaussian peak
%g = @(xi) cos(2*pi*0.5 + 0.5*xi); %Oscillatory

% %Some special G-function
%a = (1 - 2)/2.0;
%g = @(xi) (abs(4*xi - 2) + a)/(1 + a); 


%g = @(xi) abs(xi);

error = [];
error_mean = [];
error_denom = [];
polyBasis = 'Hermite';
%Gauss-Hermite/Legendre quadrature points (Q-points)
Q = 500;
[xi,w] = gaussQuad(Q,polyBasis);
g_pts = g(xi);


%MSE quadrature points 
Q_mse = 500;
[xi_mse,w_mse] = gaussQuad(Q_mse,polyBasis);


%Number of terms in the polynomial approximation
%order = [1,4,9,14,19,24,29];
order = 1:20;
%dlmwrite('coeff.dat',order);
N_iter = 1;
for N = order
    
    gamma = factorial(0:N);
    g_hat = zeros(round(N),1); 
    g_approx = 0;
    for i=1:N+1
        polynomial = hermite(xi, i-1);
        gamma_num = (sum(w.*polynomial.*polynomial));
        numer = sum(w.*polynomial.*g_pts);
        g_hat(i,N_iter) = numer/gamma(i);
        error_denom = [error_denom;abs((gamma(i) - gamma_num)/numer)];
        %g_hat(i,1) = sum(w.*polynomial.*g_pts)/(gamma(i));
        g_approx = g_approx + g_hat(i,N_iter)*hermite(xi_mse,i-1);
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

plot(xi,g_pts,'bo-','LineWidth',2);
hold on;
plot(xi_mse, g_approx,'r','LineWidth',2);
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(g_pts)-0.5,max(g_pts)+0.5]);
figure(2)
loglog(order, error, '-o','LineWidth',2,'MarkerSize',8);
hold on;
%p = polyfit(log(order),log(error)',1);
%loglog(order,order.^(-0.5*p(1)),'r','LineWidth',2);
xlabel('N');ylabel('||g(\xi) - g_N(\xi)||_{L_2}')
