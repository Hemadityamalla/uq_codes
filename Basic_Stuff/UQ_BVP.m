clear;clc;close all; format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

N = 10;
dx = 1.0/N;
x = [0:dx:1]';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
T_left = 1;
c0 = 1.0 ;
% % Define phi
d = 10; sigma = 1 ;


mc_mean = [];

n = 5*logspace(1,4,100);
c = parula(length(n));c_iter = 1;
for sim = n
    for jj=[1:1:sim]
        xi = randn(d,1);
        gmid = 0.0 ;
        for i =1:d
            gmid = gmid + xi(i)*((sqrt(2)*sigma)./((i-0.5)*pi))*(sin(pi*xmid*(i-0.5)));
        end
        cmid = c0*exp(gmid);
        c_left = cmid(2:N); % note this is because of how Matlab's spdiags works
        c_right = cmid(1:N-1);
        c_diag = -(c_left+c_right);
        A = 1/(dx^2)*spdiags([c_left c_diag c_right],[-1 0 1],N-1,N-1);
        b = zeros(N-1,1) ;
        b(1) = -(cmid(1)/(dx^2))*T_left;
        u = A\b;
        
        if ~mod(jj,5e3)
            plot(x,[1;u;0],'Color',c(c_iter,:));
            hold on;
            c_iter = c_iter+1;
        end
        T_half(jj,1) = u(N/2);
            
    end
    mc_mean = [mc_mean, mean(T_half)];
end
figure(2)
semilogx(n,mc_mean,'b.-');
hold on;
mean_temp = 0.5;
semilogx(n,mean_temp*ones(length(mc_mean),1),'r');
xlabel('n');ylabel('Mean T(x=0.5)');legend('Monte Carlo mean','Exact mean');
xpos = 500;ypos = 500; width = 800; height = 800;
set(gcf,'Position',[xpos ypos width height])
figure(3)
xpos = 500;ypos = 500; width = 800; height = 800;
set(gcf,'Position',[xpos ypos width height])
error = abs(mc_mean - mean_temp);
loglog(n,error,'ro-')
hold on;
coeffs = polyfit(log(n),log(error),1);
y = polyval(coeffs,log(n));
loglog(n,exp(y),'b');
xlabel('n'); ylabel('Absolute error');
legend('Absolute error ',strcat('line with slope ',num2str(coeffs(1))));