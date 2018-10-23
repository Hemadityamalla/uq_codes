% clear;clc;
% 
% N = 10;
% dx = 1.0/(N);
% x = [0:dx:1]; xmid = (x(1:length(x)-1) + x(2:end))/2.0;
% %Define phi 
% d = 100;sigma = 1;
% j = 1;
% for sim=[1:5:5000]
% %Define random variable coeffs
% xi = randn(d,1);
% 
% %Define g(x)
% c0 = 1.0;
% g = 0.0;
% gmid = 0.0;
% for i=1:d
%    g = g + xi(i)*((sqrt(2)*sigma)./((i-0.5)*pi))*(sin(pi*x*(i-0.5)));
%    gmid = gmid + xi(i)*((sqrt(2)*sigma)./((i-0.5)*pi))*(sin(pi*xmid*(i-0.5)));
% end
% c = c0*exp(g); %(x+1) %Test cases
% cmid = c0*exp(gmid); %(xmid+1)
% A = -2.0*diag(cmid(1:N-1)) + diag(c(2:N-1),-1) + diag(c(2:N-1),1);
% b = zeros(N-1,1);
% b(1,1) = -c(1);
% u = A\b;
% T_half(j,1) = u(5); 
% mc_mean(j,1) = mean(T_half);
% j = j+1;
% % figure(1)
% % plot(x,[1;u;0],'-');
% % xlim([0,1]);ylim([-0.5,1.5]);
% % hold on;
% end
% 
% figure(2)
% plot(mc_mean);


clear; clc;
N = 10;
dx = 1.0/N;
x = [0:dx:1]';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
T_left = 1;
c0 = 1.0 ;
% Define phi
d = 10; sigma = 1 ;
j = 1;
for sim = 1 : 5 : 1000
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
T_half(j,1) = u(N/2);
mc_mean(j,1) = mean(T_half);
j = j+1;
figure(1)
plot(x,[1;u;0]);
xlim([0,1]); ylim([-0.5,1.5]);
hold on ;
end
figure(2)
plot(mc_mean);