clear;
clc;

N = 100;
dx = 1.0/(N);
x = [0:dx:1];
xmid = (x(1:length(x)-1) + x(2:end))/2.0;
%Define phi 
d = 1000;
sigma = 1;
for sim=[1 2 3 4 5]
%Define random variable coeffs
xi = randn(d,1);

%Define g(x)
c0 = 1.0;
g = 0.0;
gmid = 0.0;
for i=1:d
   g = g + xi(i)*((sqrt(2)*sigma)./((i-0.5)*pi))*(sin(pi*x*(i-0.5)));
   gmid = gmid + xi(i)*((sqrt(2)*sigma)./((i-0.5)*pi))*(sin(pi*xmid*(i-0.5)));
end
c = c0*exp(g);
cmid = c0*exp(gmid);
A = -2.0*diag(cmid(1:N-1)) + diag(c(2:N-1),-1) + diag(c(2:N-1),1);
b = zeros(N-1,1);
b(1,1) = -c(1);

u = A\b;
plot(x,[1;u;0],'-o');
hold on;
end