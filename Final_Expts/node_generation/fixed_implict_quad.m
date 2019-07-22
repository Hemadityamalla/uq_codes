


function [x,w] = fixed_implict_quad(N, y)
    Kmax = length(y);
    %Initialize quad rule
    x = y(1:N);
    w = ones(N,1)/(N);
    %Implicit quad rule
    for D=N:(Kmax-1)
       %Node addition
       x = [x;y(D+1)];
       w = [((D)/(D+1))*w;1./(D+1)];
       %Update weights
       %V = fliplr(vander(x))';
       V = general_vandermonde(x, @(x,k) x.^(k-1), 1:N);
       nullVec = null(V(1:end,:)); 
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