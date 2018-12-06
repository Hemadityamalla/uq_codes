clear; clc;
N_points = 10;
dx = 1.0/N_points;
x = [0:dx:1]';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
T_left = 1;
c0 = 1.0 ;
sigma = 1 ;

%Gauss-Hermite quadrature points (Q points)
Q = 16;
[xi,w] = gaussQuad(Q,'Legendre');

%Deterministic solver

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

%Evaluating the coefficients and assembling the final function
ip = xi;%(xi(1):0.01:xi(end))'; %Evaluation points
error = [];
N = 2;
%for N = [1,4,9,19,24] %Number of terms in the polynomial approximation
T_hat = zeros(N+1,1); %Expansion coefficients array
T_approx = 0; %Final approximated function variable
for i=1:N+1
    polynomial = legendre(xi, i-1);
    T_hat(i,1) = sum(w.*polynomial.*T_half)/(sum(w.*polynomial.*polynomial));
    T_approx = T_approx + T_hat(i,1)*legendre(ip,i-1);
end
%error = [error;sqrt(sum(w.*(T_half - T_approx).^2))];
%end
T_mean = 0.497159169844595; %Computed using 1 million MC samples.
plot(xi,T_half,'b.','Linewidth',1,'MarkerSize',8);
hold on;
plot(ip, T_approx,'LineWidth',1);
xlim([min(xi)-0.05,max(xi)+0.05]);
ylim([min(T_half)-0.5,max(T_half)+0.5]);