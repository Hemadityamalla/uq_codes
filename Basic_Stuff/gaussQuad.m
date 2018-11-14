function [xi,w] = gaussQuad(n, type)
%Gaussian Quadrature for a given set of orthogonal polynomials
%Based on the paper of Golub and Welsch
%Generates 'n' quadrature points and the corresponding weights as two
%separate column vectors
% type = 1 corresponds to Hermite(probabilist) polynomials

if nargin == 1
    type = 'Hermite';
end

if strcmp(type,'Hermite')
    alpha = zeros(n,1);
    beta = sqrt([1:n-1])';
    L = spdiags(beta, -1, n, n);
    J = L + L' + spdiags(alpha, 0, n, n);
    [V,D] = eig(full(J));
end

xi = diag(D); %Converting the diagonal matrix into a column vector
w = V(1,:)'.^2; %Taking the first component from each eigenvector and squaring it.



end