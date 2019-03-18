clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

Kmax = 10000;

%Samples
y = rand(Kmax,1);

%Degree/ cardinality of the basis
for N = [5:5:30]
    mean_error = 0.0;
for iter = 1:10
    %QRule1
    %f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2)); coeffs = 1:N;
    %QRule2
    %f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(x.^(k/2 + 1/3)); coeffs = 1:N;
    %QRule3
    %f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(x.^(k/2 - 2/3)); coeffs = 1:N;
    %QRule4
    %f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(sin(x).*x.^(k/2)); coeffs = 1:N;
    %QRule5
    f = @(x,a) exp(x.*(a-1)); coeffs = [1,rand(1,N-1)]; coeffs = [1,rand(1,N-1)];

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

    %1) Integration with monomial of order N-1
    %exact = 1/(N);
    %approx=dotprod(x'.^(N-1),w);

    %2) Integration with cont. Genz function
    %exact = 1 - exp(-0.25);
    %approx = dotprod(0.25*exp(-0.5*abs(x' - 0.5)),w);

    %3) Integration with disc. Genz function
    %exact = 2*(exp(0.5) - exp(0.25));
    %approx = dotprod((x' < 0.5)*0 + (x' >= 0.5).*exp(0.5*x'), w);

    %4) Integration with corner peak Genz function
    exact = 1./6;
    approx = dotprod(1./(1 + 5*x').^2, w);
    
    %5) Integration with corner peak Genz function
%      exact = -1./2116;
%      approx = dotprod((x'.^45).*log(x'), w);
    
    %6) Integration with highly oscillatory function
%     exact = (2./55)*sin(55./2)^2;
%     approx = dotprod(sin(55*x'),w);
    mean_error = mean_error + abs(exact - approx);
end
     mean_error = mean_error/iter
%     if N == 5
%         dlmwrite('Errors_quad6_fn6.dat',mean_error);
%     else
%         dlmwrite('Errors_quad6_fn6.dat',mean_error,'-append');
%     end
end
 
%3) Integration with osc. Genz function
    %exact = 2*(sin(0.6*pi + 0.5) - sin(0.6*pi))
    %dotprod(cos(0.6*pi + 0.5*x'),w)
