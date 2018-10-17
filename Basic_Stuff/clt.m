clear;clc;
N = 10000; %Number of samples of the 'sampling mean'

%Generate exponentially distributed random numbers
lambda = 10; %Mean and standard deviation is 1/lambda
mu = 1/lambda; sigma = 1/lambda;
n = 100; %number of samples
for i=1:N
   Xn(i,1) = mean(exprnd(mu, [n,1])); 
end
%Computing the mean and std. dev. of Xn
mu_Xn = mu;
sigma_Xn = sigma/sqrt(n);
%Plotting 
histogram(Xn,'Normalization','pdf');
hold on;
x = linspace(mu_Xn - 5*sigma_Xn,mu_Xn + 5*sigma_Xn,5000);
plot(x,(1/sqrt(2*pi*sigma_Xn^2))*exp(-(x - mu).^2./(2*sigma_Xn^2)),'r');