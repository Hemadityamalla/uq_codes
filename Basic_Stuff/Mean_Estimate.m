clear;
clc;

%Estimating the mean of a random variable

for i=10.^[2 3 4 5 6 7]
Y = randn(i,1);
X = exp(-1+2*Y);
mean(X)
end