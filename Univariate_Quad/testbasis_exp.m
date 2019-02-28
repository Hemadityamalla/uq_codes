clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

N = 5;
a= N;
f = @(x,k) exp(a*x).*(x).^k;
% g = @(x,k) exp(-a*x).*(x).^k;

%Generate the set of points
x = cos((pi*(0:N-1))/(N-1)); %cos((pi*(k-1))/(N-1))

%Defining the weighting function 
%w = @(x) 0.5; %Uniform distribution

%Vandermonde matrix
V = vander_analogue(N, x, f);
%V = vander_cheb(N,x');


%RHS integral vector
Eint = @(x,k) (exp(a*x)/a^(k+1)).*(((-1).^(0:k)).*(factorial(0:k))*((a*x).^(k-(0:k)))');
b = zeros(N,1);
b(1,1) = N; 
for idx = 2:N
    b(idx,1) = (Eint(1,idx-1)-Eint(-1,idx-1));
end

%Weights
w = V\b;