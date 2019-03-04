clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

n=4;

f = @(x,a) exp(a*x); %Here 'a' varies and the basis set has 'a' that are pairwise distinct

%Generate the set of points
x = -1+2*rand(1,n);%cos((pi*(0:n-1))/(n-1));

%Defining the weighting function 
%w = @(x) 0.5; %Uniform distribution


coeff = 0:(2*n-1);

%Checking if the given system is complete-Chebyshev on [-1,1]
det(f(x, (0:length(x)-1)'))

%Vandermonde matrix
V = zeros(2*n, n);
for ii=1:(2*n)
   V(ii,:) = f(x,coeff(ii)); 
end

%RHS integral vector
b = zeros(2*n,1);
b(1,1) = 2.0; %The first basis function is 1 
for idx = 2:(2*n)
    b(idx,1) = (exp(coeff(idx)) - exp(-coeff(idx)))/coeff(idx); 
end

%Weights
w = V\b