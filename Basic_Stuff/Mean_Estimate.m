clear;clc;close all; format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

%Estimating the mean of a random variable y = f(X)
j=1;
n = logspace(1,6);
for i=n
    Y = rand(round(i),1);
    X = exp(-1+2*Y);
    mu_y(j,1) = mean(X);
    j = j+1;
end
exact_mean = (exp(1) - exp(-1))/2;
figure(1);
plot(n,mu_y,'b.-');
hold on;
plot(exact_mean*ones(1e6,1),'r');
legend('estimated mean','exact mean');
xlabel('n (Number of samples)'); ylabel('Mean of Y');
figure(2);
error = abs(mu_y - exact_mean);
loglog(n,error,'ro-');
hold on;
coeffs = polyfit(log(n),log(error'),1);
y = polyval(coeffs,log(n));
loglog(n,exp(y),'b');
xlabel('n (Number of samples)'); ylabel('Absolute error');
legend('Absolute error ',strcat('line with slope',num2str(coeffs(1))));