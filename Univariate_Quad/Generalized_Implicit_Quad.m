clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

Kmax = 10000;

%Samples
y = rand(Kmax,1);

%Degree/ cardinality of the basis
N = 10;
%f = @(x,a) x.^(0.5*(a-1)); %Here 'a' varies and the basis set has 'a' that are pairwise distinct
%f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(x.^(k/2 + 1/3)); %this must have the coeffs starting from 1
%f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2));
%f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(sin(x).*x.^(k/2));
f = @(x,a) exp(x.*(a-1));
coeffs = [1,rand(1,N-1)];

%Initialize quad rule
x = y(1:N);
w = ones(N,1)/(N);

%Implicit quad rule
for D=N:(Kmax-1)
   %Node addition
   x = [x;y(D+1)];
   w = [((D)/(D+1))*w;1./(D+1)];
   %Update weights
   V = general_vandermonde(x, f, coeffs);
   nullVec = null(V);
   c = nullVec(:,1);
   
   % Apply Caratheodory's and determine the node that will be removed
    [alpha, k] = min(w(c>0)./c(c>0));
    %[alpha, k] = max(w(c<0)./c(c<0));
    id = find(c > 0);%find(c < 0)
    k = id(k);

    %Node removal
    w = w-alpha*c;
    x(k, :) = []; w(k) = [];
    
end

%Numerical tests

%Integration with monomial of order N-1
exact = 1/(N-1)
dotprod(x'.^(N-1),w)

%Integration with cont. Genz function