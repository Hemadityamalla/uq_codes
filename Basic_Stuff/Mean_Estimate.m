clear;clc;
%Estimating the mean of a random variable y = f(X)
j=1;
n = [1:20:2e4];
for i=n
    Y = rand(i,1);
    X = exp(-1+2*Y);
    mu_y(j,1) = mean(X);
    j = j+1;
end
exact_mean = (exp(1) - exp(-1))/2;
figure(1);
plot(mu_y,'b.-');
hold on;
plot(exact_mean*ones(length(mu_y),1),'r','LineWidth',2);
legend('estimated mean','exact mean');
figure(2);
error = abs(mu_y - exact_mean);
loglog(n,error,'ro-');
hold on;
y = polyval(polyfit(log(n),log(error'),1),log(n));
loglog(n,exp(y),'b','LineWidth',2);
xlabel('n'); ylabel('|estimated_{\mu_y} - exact_{\mu_y}|');
legend('|estimated_{\mu_y} - exact_{\mu_y}|','line with slope -1/2');