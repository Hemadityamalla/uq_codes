clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

N = 50;
%Generate the set of points
k=1:N;
x = cos((pi*(k-1))/(N-1));

%Defining the weighting function 
%w = @(x) 0.5; %Uniform distribution

%Vandermonde matrix
V = fliplr(vander(x));
%V = vander_cheb(N,x');


%RHS integral vector
b = @(x) (1./k).*(x.^k);
% b = zeros(N,1);
% [xgauss, wgauss] = gaussQuad(100,'Legendre');
% for idx = 1:N
%     %b(idx,1) =  dotprod(chebyshev(xgauss,idx-1)',wgauss);
%     b(idx,1) = 0.5*dotprod(xgauss'.^(idx-1),wgauss);
%end

%Weights
w = V\(b(1) - b(-1))';
%w = V\b;