function [xi,w] = gaussQuad(n, type)
%Gaussian Quadrature for a given set of orthogonal polynomials
%Based on the paper of Golub and Welsch
%Generates 'n' quadrature points and the corresponding weights as two
%separate column vectors
% type = 'Hermite' corresponds to Hermite(probabilist) polynomials- default

if nargin == 1
    type = 'Hermite';
end

if strcmp(type,'Hermite') %Normal distn.
    alpha = zeros(n,1);
    beta = sqrt(1:n-1)';
    L = spdiags(beta, -1, n, n);
    J = L + L' + spdiags(alpha, 0, n, n);
elseif strcmp(type,'Legendre') %Uniform distn. in [-1,1]
    alpha = zeros(n,1);
    beta = (sqrt((1:n-1).^2./(4*(1:n-1).^2 - 1)))';
    L = spdiags(beta, -1, n, n);
    J = L + L' + spdiags(alpha, 0, n, n);
    [V,D] = eig(full(J));
    %Correct this part!
    xi = diag(D); %Converting the diagonal matrix into a column vector
    w = 2*V(1,:)'.^2; %Taking the first component from each eigenvector and squaring it.
    return
elseif strcmp(type,'Chebyshev') %Chebyshev distn. in [-1,1]
    alpha = zeros(n,1);
    beta = 0.5*ones(n-1,1);
    L = spdiags(beta, -1, n, n);
    J = L + L' + spdiags(alpha, 0, n, n);
end


[V,D] = eig(full(J));
xi = diag(D); %Converting the diagonal matrix into a column vector
w = V(1,:)'.^2; %Taking the first component from each eigenvector and squaring it.



end