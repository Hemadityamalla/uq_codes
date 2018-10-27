clear;
clc;

N = 10000;
a = 1;
b = 1;
mu = 0;
sigma = 1;
%4 point gaussian quadrature
w = [0.652145, 0.652145, 0.347855, 0.347855];
t = [0.339981, -0.339981, 0.861136, -0.861136];


X = mu + sigma*randn(N,1);
f(1,1) = 0.5*dot(w,betainv((t+1)/2.0,a,b)*1.0);
f(2,1) = 0.5*dot(w,betainv((t+1)/2.0,a,b).*norminv((t+1)/2.0));
f(3,1) = 0.5*dot(w,betainv((t+1)/2.0,a,b).*(norminv((t+1)/2.0).^2 - 1));
f(4,1) = 0.5*dot(w,betainv((t+1)/2.0,a,b).*(norminv((t+1)/2.0).^3 - 3*norminv((t+1)/2.0)));
f(5,1) = 0.5*dot(w,betainv((t+1)/2.0,a,b).*(norminv((t+1)/2.0).^4 - 6*norminv((t+1)/2.0).^2 + 3));
f(6,1) = 0.5*dot(w,betainv((t+1)/2.0,a,b).*(norminv((t+1)/2.0).^5 - 10*norminv((t+1)/2.0).^3 + 15*norminv((t+1)/2.0)));

Y_n = f(1,1)+ f(2,1)*X + 0.5*f(3,1)*(X.^2 - 1) + f(4,1)/6.0*(X.^3 - 3*X) + f(5,1)/24.0*(X.^4 - 6*X.^2 + 3) + f(6,1)/120.0*(X.^5 - 10*X.^3 + 15*X);

[p,x] = ksdensity(Y_n);plot(x,p);
hold on;
[q,x] = ksdensity(betarnd(a,b,[N,1]));plot(x,q);
legend('approximate','exact')
