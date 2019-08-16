
%f_deets.fn = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2))
%f_deets.coeffs = 1:N;
%[x,w] = general_fixed_implicit_quad(f_deets, N, Kmax);


function [x,w,y] = general_fixed_implicit_quad(f_deets, N, Kmax)
    %Initialize the coefficients of the basis fucntions (if not initialized)
    if isempty(f_deets.coeffs)
        f_deets.coeffs = 1:N;
    end
    %Initialize quad rule
    y = rand(Kmax,1);%gamrnd(2,0.5,[Kmax,1]);
    x = y(1:N);
    w = ones(N,1)/(N);
    %Implicit quad rule
    for D=N:(Kmax-1)
       %Node addition
       x = [x;y(D+1)];
       w = [((D)/(D+1))*w;1./(D+1)];
       %Update weights
       V = general_vandermonde(x, f_deets.fn, f_deets.coeffs);
       nullVec = null(V);
       c = nullVec(:,1);

       % Apply Caratheodory's and determine the node that will be removed
        [alpha, k] = min(w(c>0)./c(c>0));
        id = find(c > 0);
        if isempty(alpha)
            [alpha, k] = max(w(c<0)./c(c<0));
            id = find(c < 0);
        end
        
        k = id(k);

        %Node removal
        w = w-alpha*c;
        x(k, :) = []; w(k) = [];
    end

end