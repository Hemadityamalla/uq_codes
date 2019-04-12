clear;clc; format long;

%N-- desired degree of the quadrature

N = 5;

%Sample Integrand
u = @(x) (x <= 0.5)*0 + (x > 0.5).*exp(0.5*x);
exact = -2*exp(-0.5) + 2*exp(-0.25);

%Generating random samples
x = rand(N,1);

%Initial quadrature rule
xquad = 0.5; wquad = 1; %This is based on clenshaw-curtis
V = zeros(N,N);
for ii=1:4
   %Evaluation
   fn = u(xquad);
   if length(xquad) == 1
       V(1,:) = ones(N,1)*fn;
   end
   %Adding the new basis
   V(ii,:) = 
end

