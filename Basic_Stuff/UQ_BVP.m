clear; clc;
N = 100;
dx = 1.0/N;
x = [0:dx:1]';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
T_left = 1;
c0 = 1.0 ;
% % Define phi
d = 10; sigma = 1 ;


j = 1;

%One rv diffusion constant- log normal
mu = 0;
sigma = 0.75;

n = [1:5:10e3];
for sim = n
xi = randn(d,1);%randn(1,1)*ones(length(xmid),1);
gmid = 0.0 ;
for i =1:d
gmid = gmid + xi(i)*((sqrt(2)*sigma)./((i-0.5)*pi))*(sin(pi*xmid*(i-0.5)));
end
cmid = exp(gmid);
c_left = cmid(2:N); % note this is because of how Matlab's spdiags works
c_right = cmid(1:N-1);
c_diag = -(c_left+c_right);
A = 1/(dx^2)*spdiags([c_left c_diag c_right],[-1 0 1],N-1,N-1);
b = zeros(N-1,1) ;
b(1) = -(cmid(1)/(dx^2))*T_left;
u = A\b;
%c(j,1) = cmid(N/2);
T_half(j,1) = u(N/2);
mc_mean(j,1) = mean(T_half);
j = j+1;
 %figure(1)
 %plot(cmid,'LineWidth',2);
% xlim([0,1]); ylim([-0.5,1.5]);
 %hold on ;
end
figure(2)
plot(mc_mean,'b.-','LineWidth',1);
hold on;
mean_temp = 0.49999;
plot(mean_temp*ones(length(mc_mean),1),'r','LineWidth',2);
figure(3)
error = abs(mc_mean - mean_temp);
loglog(n,error,'ro-')
hold on;
coeffs = polyfit(log(n),log(error'),1);
y = polyval(coeffs,log(n));
loglog(n,exp(y),'b','LineWidth',2);
xlabel('n'); ylabel('|estimated_{\mu_y} - exact_{\mu_y}|');
legend('|estimated_{\mu_y} - exact_{\mu_y}|',strcat('line with slope',num2str(coeffs(1))));