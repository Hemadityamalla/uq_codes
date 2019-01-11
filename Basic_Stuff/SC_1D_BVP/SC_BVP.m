clear;clc;close all; format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

N_points = 12;
dx = 1.0/N_points;
x = [0:dx:1]';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
T_left = 1;
c0 = 1.0;
sigma = 1;
polyBasis = 'Hermite';
%Gauss-Hermite quadrature points (Q points)
Q = 100;
[xi,w] = gaussQuad(Q,polyBasis);

%Deterministic solver
% Q = 1e7;
% xi = -1 + 2*rand(Q,1);
for i=[1:1:Q]
    gmid = xi(i)*((sqrt(2)*sigma)./((0.5)*pi))*(sin(pi*xmid*(0.5)));
    cmid = c0*exp(gmid);
    c_left = cmid(2:N_points); % note this is because of how Matlab's spdiags works
    c_right = cmid(1:N_points-1);
    c_diag = -(c_left+c_right);
    A = 1/(dx^2)*spdiags([c_left c_diag c_right],[-1 0 1],N_points-1,N_points-1);
    b = zeros(N_points-1,1) ;
    b(1) = -(cmid(1)/(dx^2))*T_left;
    u = A\b;
    T_half(i,1) = u(N_points/2);
end
T_mean_normal = 0.497138364518936;
T_mean_uniform = 0.498927930331782;
T_1_Normal_Hermite = [0.5;0.496754747025497;0.497076796853061;0.497080433744812;0.497080434486761];
T_1_Normal_Legendre = [0.5;0.498874098722895;0.498891853597968;0.498891854987300;0.498891854987301];
T_1_Uniform_Hermite = [0.5;0.496754747025497;0.497076796853061;0.497080433744812;0.497080434486761];
T_1_Uniform_Legendre = [0.5;0.498874098722895;0.498891853597968;0.498891854987300;0.498891854987301];

%Evaluating the coefficients and assembling the final function
error = [];
order = 1:29;
for N = order %Number of terms in the polynomial approximation
T_hat = zeros(N+1,1); %Expansion coefficients array
T_approx = 0; %Final approximated function variable
for i=1:N+1
    %Pre-computing the normalization factors
    if strcmp(polyBasis,'Hermite')
        gamma = factorial(0:N); %Hermite
    else
        gamma = 2.0./(2*(0:N) + 1.0); %Legendre
    end
    polynomial = hermite(xi, i-1);
    T_hat(i,1) = sum(w.*polynomial.*T_half)/(gamma(i));
    T_approx = T_approx + T_hat(i,1)*hermite(xi,i-1);
    
end
error = [error;sqrt(sum(w.*(T_half - T_approx).^2))];
    hold on;
    plot(xi, T_approx,'LineWidth',2);
    xlim([min(xi)-0.05,max(xi)+0.05]);
end
figure(1)
plot(xi,T_half,'b.','Linewidth',1,'MarkerSize',8);
hold on;
plot(xi, T_approx,'LineWidth',1);
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(T_half)-0.5,max(T_half)+0.5]);
figure(2)
loglog(order, error, '-o');
grid on
hold on;
%p = polyfit(log(order),log(error)',1);
%loglog(order,order.^(-0.5*p(1)),'r');
xlabel('N (Number of gPC expansion terms)');ylabel('Mean-square Error');