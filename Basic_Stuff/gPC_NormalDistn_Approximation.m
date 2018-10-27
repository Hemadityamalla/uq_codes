clear;
clc;

N = 5000;
mu = 0;
sigma = 1;

%4 point gaussian quadrature
w = [0.652145, 0.652145, 0.347855, 0.347855];
t = [0.339981, -0.339981, 0.861136, -0.861136];


Z = randn(N,1);


FyInv = norminv((t+1)/2.0,mu,sigma);
FxInv = norminv((t+1)/2.0);
f(1,1) = 0.5*dot(w,FyInv*1.0);
f(2,1) = 0.5*dot(w,FyInv.*FxInv);
f(3,1) = 0.5*dot(w,FyInv.*(FxInv.^2 - 1));
f(4,1) = 0.5*dot(w,FyInv.*(FxInv.^3 - 3*FxInv));
f(5,1) = 0.5*dot(w,FyInv.*(FxInv.^4 - 6*FxInv.^2 + 3));
f(6,1) = 0.5*dot(w,FyInv.*(FxInv.^5 - 10*FxInv.^3 + 15*FxInv));

Y_n = f(1,1)+ f(2,1)*Z + f(3,1)/2.0*(Z.^2 - 1) + f(4,1)/6.0*(Z.^3 - 3*Z) + f(5,1)/24.0*(Z.^4 - 6*Z.^2 + 3) + f(6,1)/120.0*(Z.^5 - 10*Z.^3 + 15*Z);

%[p,x] = ksdensity(Y_n);plot(x,p);
histogram(Y_n,'Normalization','pdf');
hold on;
histogram(mu + sigma*randn(N,1),'Normalization','pdf');
%[q,x] = ksdensity(mu + sigma*randn(N,1));plot(x,q);
legend('approximate','exact')
