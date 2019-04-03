%clear;format long;
%set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

function main(Qrule, testFn)

Kmax = 1000;
% Qrule = 1;
% testFn = 2;
switch testFn
    case 1
        u = @(x) x.^(9);
    case 2
        u = @(x) 0.25*exp(-0.5*abs(x - 0.5));
    case 3
        u = @(x) 1./(1 + 5*x).^2;
    case 4
        u = @(x) x.^(45).*log(x);
    case 5
        u = @(x) sin(55*x);
end

%Degree/ cardinality of the basis
for N = [5:5:30]
    mean_error = 0.0;
    switch Qrule
    case 1
    %QRule1
    f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2)); coeffs = 1:N;
    case 2
    %QRule2
    f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(x.^(k/2 + 1/3)); coeffs = 1:N;
    case 3
    %QRule3
    f = @(x,a) exp(x.*(a-1)); coeffs = [1,rand(1,N-1)]; coeffs = [1,rand(1,N-1)];
    case 4
    %QRule4
    f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(sin(x).*x.^(k/2)); coeffs = 1:N;
    case 5
    %QRule5
    f = @(x,k) 1./(1 + k*x); coeffs = 1 - 1./sqrt(1:N);
    %f = @(x,k) x.^(k-1); coeffs = 1:N;
    end
Navg = 50;
for iter = 1:Navg
    %Samples
    if  Qrule == 5
        y = -1 + 2*rand(Kmax,1);
    else
        y = rand(Kmax,1);
    end
    %Initialize quad rule
    x = y(1:N);
    w = ones(N,1)/(N);

    %Implicit quad rule
    for D=N:(Kmax-1)
       %Node addition
       x = [x;y(D+1)];
       w = [((D)/(D+1))*w;1./(D+1)];%Update weights
       V = general_vandermonde(x, f, coeffs);
       %nullVec = null(V);
       %size(nullVec)
       [~,~,U] = svd(V);
       nnv = size(V,2) - size(V,1);
       nullVec = U(:,end-nnv:end);
       c = nullVec(:,end);

       % Apply Caratheodory's and determine the node that will be removed
        [alpha, k] = min(w(c>0)./c(c>0));
        %[alpha, k] = max(w(c<0)./c(c<0));
        id = find(c > 0);%find(c < 0)
        k = id(k);

        %Node removal
        w = w-alpha*c;
        x(k, :) = []; w(k) = [];

    end


    exact = mean(u(y'));
    approx = dotprod(u(x'),w);
    mean_error = mean_error + abs(exact - approx);
end
     mean_error = mean_error/iter;
    fname = strcat('Errors_quad',num2str(Qrule),'_fn',num2str(testFn),'.dat');
    if N == 5
        dlmwrite(fname,mean_error);
    else
        dlmwrite(fname,mean_error,'-append');
    end
end

end