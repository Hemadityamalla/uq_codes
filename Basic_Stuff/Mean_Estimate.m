clear;clc;
%Estimating the mean of a random variable y = f(X)
j=1;
for i=[1:20:2e4]
Y = rand(i,1);
X = exp(-1+2*Y);
mu_y(j,1) = mean(X);
j = j+1;
end
plot(mu_y,'b.-');
hold on;
plot(0.5*(exp(1)-exp(-1))*ones(length(mu_y),1),'r','LineWidth',2);