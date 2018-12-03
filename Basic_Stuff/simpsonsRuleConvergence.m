clear; clc;
%Script to check the order of covergence of the composite trapezoidal rule

f = @(x) x;

a = 0;
b = 1;

n = 10*[1,2,4,8,16,32,64];
In = [];
for i=n
    x = linspace(a,b,i);
    In = [In,((b-a)/i)*(sum(f(x))-0.5*(f(a) + f(b)))];
end

I = 0.5;

error = abs(I - In);
h = (b-a)./n;

polyfit(log(h),log(error),1)