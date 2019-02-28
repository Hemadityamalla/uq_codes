clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

N = 50;
%Generate the set of points
k=1:N;
x = cos((pi*(k-1))/(N-1));

%Vandermonde matrix
V = fliplr(vander(x))';
%V = vander_cheb(N,x');


%RHS integral vector ( (x^n)/(n+1) )

b = zeros(N,1);
b(:,1) = (1./k).*((1).^(k) - (-1).^(k));


%Weights
w = V\b;