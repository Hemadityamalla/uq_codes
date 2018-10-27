clear;
clc;

N = 5000;

mu = 0;
sigma = 1;
Z = randn(N,1);
X = mu + sigma*Z;

histogram(X,'Normalization','pdf');

Y = exp(X); %Exact function 
Yn = exp(mu + 0.5*sigma^2)*(1 + sigma*Z + 0.5*sigma^2*(Z.^2 - 1) + (sigma^3/3)*(Z.^3 - 3*Z)); %Approximate function
nbins = 50;
histogram(Yn,nbins,'Normalization','pdf');

y = linspace(0,20,N);
hold on;
%The pdf below is usually known
plot(y,1.0./(y*sigma*sqrt(2*pi)).*exp(-(log(y) - mu).^2/(2*sigma^2)),'LineWidth',2);
