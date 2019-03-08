clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

m = 500;
[x,w] = clencurt(m);
x = 0:(1/(m)):1;%(x + 1)*0.5;
%Defining the basis function

%f = @(x,k) (mod(k,2)==0)*(x.^(k/2)) + (mod(k,2)==1)*(log(x).*x.^((k-1)/2));
%f = @(x,k) (mod(k-1,2)==0)*(x.^((k-1)/2)) + (mod(k-1,2)==1)*(x.^(1/3).*x.^((k)/2));
f = @(x,k) (mod(k-1,2)==0)*(x.^((k-1)/2)) + (mod(k-1,2)==1)*(x.^(1/3).*x.^((k)/2));


V = general_vandermonde(x,f,1:m+1);