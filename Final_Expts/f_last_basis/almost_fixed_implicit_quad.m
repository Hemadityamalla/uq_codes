%This is trivial and wrong as we are just adding the function evaluations
%and trying to obtain a quadrature rule for it.

function [x,w,y] = almost_fixed_implict_quad(N,f, Kmax)
    y = rand(Kmax,1); %Uniform distribution
    %Initialize quad rule
    x = y(1:N);
    w = ones(N,1)/(N);
    %Implicit quad rule
    for D=N:(Kmax-1)
       %Node addition
       x = [x;y(D+1)];
       w = [((D)/(D+1))*w;1./(D+1)];
       %Update weights
       V = general_vandermonde(x, @(x,k) x.^(k-1), 1:length(x));
       %Replacing the second row with an interpolation of the function
       V(end,:) = f(x); %This is a bad idea as x keeps changing every iteration and each function evaluation is really expensive

       nullVec = null(V); %Wait, why do I neglect the last column?!
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

end