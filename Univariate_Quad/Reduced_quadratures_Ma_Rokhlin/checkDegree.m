function n = checkDegree(x,w,integrand,integral)
if nargin == 2
    integral = @(k) 1.0/k; %A default polynomial case
    integrand = @(x,k) x.^(k-1);
end


N = length(w); 
accuracy = 1e-14; %This causes the output to be +/- 1 or 2 of the actual degree
%Checking for degree upto 2*N
for n = 1:(2*N)
    %abs(dotprod(integrand(x',n),w) - integral(n))
    if abs(dotprod(integrand(x',n),w) - integral(n)) > accuracy
        
        return 
    end
    
end



end