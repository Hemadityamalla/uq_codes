%clear; clc;
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

%Gauss-Hermite quadrature points
%(2 point)
%xi = sqrt(2)*[-1/sqrt(2), 1/sqrt(2)];
%w = [0.88622693, 0.88622693];
%(3 point)
%xi = sqrt(2)*[-1.22474487,0.,1.22474487];
%w = [0.29540898, 1.1816359 , 0.29540898];
%(10 point)
%xi = sqrt(2)*[-3.43615912, -2.53273167, -1.75668365, -1.03661083, -0.34290133,0.34290133,  1.03661083,  1.75668365,  2.53273167,  3.43615912];
%w = [7.64043286e-06, 1.34364575e-03, 3.38743945e-02, 2.40138611e-01,6.10862634e-01, 6.10862634e-01, 2.40138611e-01, 3.38743945e-02,1.34364575e-03, 7.64043286e-06];
%(20 point)
xi = sqrt(2)*[-5.38748089, -4.60368245, -3.94476404, -3.34785457, -2.78880606,...
                -2.254974  , -1.73853771, -1.23407622, -0.73747373, -0.24534071,...
                0.24534071,  0.73747373,  1.23407622,  1.73853771,  2.254974,...
                2.78880606,  3.34785457,  3.94476404,  4.60368245,  5.38748089];
w = [2.22939365e-13, 4.39934099e-10, 1.08606937e-07, 7.80255648e-06,...
       2.28338636e-04, 3.24377334e-03, 2.48105209e-02, 1.09017206e-01,...
       2.86675505e-01, 4.62243670e-01, 4.62243670e-01, 2.86675505e-01,...
       1.09017206e-01, 2.48105209e-02, 3.24377334e-03, 2.28338636e-04,...
       7.80255648e-06, 1.08606937e-07, 4.39934099e-10, 2.22939365e-13];
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

%Hermite polynomials(probabilist)
h0 = @(xi) ones(size(xi));
h1 = @(xi) xi;
h2 = @(xi) xi.^2 - 1.0;
h3 = @(xi) xi.^3 - 3.0*xi;
h4 = @(xi) xi.^4 - 6.0*xi.^2 + 3.0;


T_0 = sum(w.*h0(xi).*T_half')/(factorial(0)*sqrt(pi));
T_1 = sum(w.*h1(xi).*T_half')/(factorial(1)*sqrt(pi));
T_2 = sum(w.*h2(xi).*T_half')/(factorial(2)*sqrt(pi));
T_3 = sum(w.*h3(xi).*T_half')/(factorial(3)*sqrt(pi));


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