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

%Gaussian quadrature points(5 point)
xi = [-2.020184, -0.95857, 0, 0.95857, 2.020184];
w = [1.18148, 0.98658, 0.94530, 0.98658, 1.18148];
Q = length(xi);

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

%Hermite polynomials
h0 = @(xi) ones(size(xi));
h1 = @(xi) 2*xi;
h2 = @(xi) 4*xi.^2 - 2.0;
h3 = @(xi) 8*xi.^3 - 12.0*xi;


T_0 = sum(gaussmf(xi,[1 0]).*w.*h0(xi).*T_half')/factorial(0);
T_1 = sum(gaussmf(xi,[1 0]).*w.*h1(xi).*T_half')/factorial(1);
T_2 = sum(gaussmf(xi,[1 0]).*w.*h2(xi).*T_half')/factorial(2);
T_3 = sum(gaussmf(xi,[1 0]).*w.*h3(xi).*T_half')/factorial(3);


% figure(2)
% plot(mc_mean,'b.-','LineWidth',1);
% hold on;
% mean_temp = 0.5;
% plot(mean_temp*ones(length(mc_mean),1),'r','LineWidth',2);
% figure(3)
% error = abs(mc_mean - mean_temp);
% loglog(n,error,'ro-')
% hold on;
% coeffs = polyfit(log(n),log(error),1);
% y = polyval(coeffs,log(n));
% loglog(n,exp(y),'b','LineWidth',2);
% xlabel('n'); ylabel('|estimated_{\mu_y} - exact_{\mu_y}|');
% legend('|estimated_{\mu_y} - exact_{\mu_y}|',strcat('line with slope',num2str(coeffs(1))));