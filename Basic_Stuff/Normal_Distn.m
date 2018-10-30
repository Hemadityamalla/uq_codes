clear;
clc;

%Generate samples that are normally distributed
N = 10000; %Number of samples
mu = 1.0; %mean
sigmasq = 4.0; %variance
X = mu + sqrt(sigmasq)*randn(N,1);

nbins = 25;
histogram(X,nbins,'Normalization','pdf');
hold on;
x = linspace(-10,10,N);
plot(x,(1/sqrt(2*pi*sigmasq))*exp(-(x - mu).^2./(2*sigmasq)));

