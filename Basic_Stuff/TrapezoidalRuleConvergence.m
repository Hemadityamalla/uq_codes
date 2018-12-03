clear; clc;
%Script to check the order of covergence of the composite trapezoidal rule

f = @(x) x.^3;%cos(x);

a = 0;
b = 1;

n = 10*[128,64,32,16,8,4,1];
In = [];
for i=n
    x = linspace(a,b,i+1);
    In = [In,((b-a)/i)*(sum(f(x))-0.5*(f(a) + f(b)))];
end

I = 0.25;

error = abs(I - In);
h = (b-a)./n;

polyfit(log(h),log(error),1) %This will give a slope approximately 2