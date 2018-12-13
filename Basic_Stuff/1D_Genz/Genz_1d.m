clear;clc;

%Defining the 1D Genz functions

g = @(xi) exp(-0.5*abs(xi)); %continuous 
%g = @(xi) 1.0./(0.5^(-2) + (xi).^2); %Product Peak 
%g = @(xi) (1 + 5*xi).^(-2); %Corner peak 
%g = @(xi) (xi > 0)*0 + (xi <= 0).*exp(0.5*xi); %Discontinuous
%g = @(xi) exp(-(0.5^2)*(xi).^2); %Gaussian peak
%g = @(xi) cos(2*pi*0.5 + 0.5*xi); %Oscillatory

% %Some special G-functionlegendre
%a = (1 - 2)/2.0;
%g = @(xi) (abs(4*xi - 2) + a)/(1 + a); 


%g = @(xi) abs(xi);

error = [];
%Gauss-Hermite/Legendre quadrature points (Q-points)
Q = 1000;
[xi,w] = gaussQuad(Q,'Hermite');
g_pts = g(xi);

ip = xi;%(xi(1):0.01:xi(end))'; %Interpolating points

%Number of terms in the polynomial approximation
order = [1,4,9,14,19,24];

for N = order
    g_hat = zeros(round(N),1); 
    g_approx = 0;
    for i=1:N+1
        polynomial = hermite(xi, i-1);
        g_hat(i,1) = sum(w.*polynomial.*g_pts)/(sum(w.*polynomial.*polynomial));
        g_approx = g_approx + g_hat(i,1)*hermite(ip,i-1);
    end
    %error = [error; abs(g_hat(1,1) - gmean)];
   %error = [error; sqrt(sum((exp(-0.5*ip.^2)/sqrt(2*pi)).*abs(g(ip) - g_approx).^2))];
   error = [error;sqrt(sum(w.*(g(ip) - g_approx).^2))];
   plot(xi,g_pts,'bo-','LineWidth',2);
hold on;
plot(ip, g_approx,'LineWidth',2);
xlim([min(xi)-0.05,max(xi)+0.05]);
end

plot(xi,g_pts,'bo-','LineWidth',2);
hold on;
plot(ip, g_approx,'r','LineWidth',2);
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(g_pts)-0.5,max(g_pts)+0.5]);
figure(2)
loglog(order, error, '-o','LineWidth',2,'MarkerSize',8);
hold on;
p = polyfit(log(order),log(error)',1);
loglog(order,order.^(-0.5*p(1)),'r','LineWidth',2);
xlabel('N');ylabel('||g(\xi) - g_N(\xi)||_{L_2}')
