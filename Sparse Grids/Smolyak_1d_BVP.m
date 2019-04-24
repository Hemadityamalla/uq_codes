clear;clc; 
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
Npts = 20;
dx = 1.0/Npts;
x = (0:dx:1)';
xmid = ( x(1:end-1) + x (2:end) ) / 2.0 ;
u_left = 1;

d = 2; %dimension of the random vector
l=(1:d)';
sigma = 1.0;
% Input for type of orthogonal polynomial basis. Alternatively, use 'Hermite'
polyBasis = 'Legendre';
quadrature = 'ClenshawCurtis';
growth = @(x) 2^(x-1);
%mu = [];1
%var = [];
% skew = [];
% kurt = [];
numpts = [];
%exact_mean = 0.529175686352757; %d=2
exact_mean = 0.586605976838356; %d=5
error_mean = [];
% exact_variance = 0.080571966133941;1
% exact_skewness = 0.285214735020665;
% exact_kurtosis = 1.493438201045323;
u_sol = [];
%[xi_mse,w_mse] = smolyakSparseGrid(d,7,growth, quadrature);
xi = -1+2*rand(d,1e6);
for k=[2,3,4,5,6,7,8]
    umean = 0;
    [xi,w] = smolyakSparseGrid(d,k,growth, quadrature);
    %Rescaling to the domain [0,1]---- not sure why I am dividing by 2^d- figure it out dumbass!.
    %xi = 0.5*(xi+1);
    %w = w/2^d;
    %u_mid = zeros(size(w));
    u_mid = [];
    for ii=1:length(xi)
    %Solution of the ODE
        %1) Computing the random diffusivity
        phi = sin((l-0.5)*pi*xmid');
        a = exp(((sqrt(2)*sigma./(l-0.5)*pi).*xi(:,ii)')*phi);
        c_left = a(2:Npts)'; % note this is because of how Matlab's spdiags works
        c_right = a(1:Npts-1)';
        c_diag = -(c_left+c_right);
        A = 1/(dx^2)*spdiags([c_left c_diag c_right],[-1 0 1],Npts-1,Npts-1);
        b = zeros(Npts-1,1) ;
        b(1) = -(a(1)/(dx^2))*u_left;
        u = A\b;
        %u_mid(1,ii) = u(Npts/2);
        u_mid = [u_mid; u(Npts/2)];
        error_mean = [error_mean; abs(mean(u_mid) - exact_mean)];
        %u_sol = [u_sol, [1;u;0]];
        %plot(x,[1;u;0]);
        %hold on;
    end
    %plot(x,umean);
    %hold on;
    degree = 2; %maximum degree of the multivariate polynomial---------experiment on this
    lexOrdering = monomialDegrees(d,degree);
    %Pre-computing the normalization factors
%     gamma = 2.0./(2*(0:degree) + 1.0); %Legendre
%     P = size(lexOrdering,1);
%     fhat = zeros(P,1);
%     fapprox = 0;
%     for i_P = 1:P
%         fhat(i_P,1) = (sum(legendre(xi',lexOrdering(i_P,:)').*u_mid.*w'))/prod(gamma(lexOrdering(i_P,:)+1));
%         %fapprox = fapprox + fhat(i_P,1)*legendre(xi_sample', lexOrdering(i_P,:)');
%     end
    fhat = gpc_coeffs(lexOrdering,xi', w'/2^d, u_mid, 'legendre'); %Scaling the weight to match the weighting function of Uniform distribution
    fapprox = gpc_polyval(fhat, xi', 'legendre', lexOrdering);
    numpts = [numpts;length(w)];
    %w = ones(1,2e6)*5e-7;
    %ap_var = dotprod((fapprox' - fhat(1,1)).^2,w');
    %ap_skew = dotprod(((fapprox' - fhat(1,1))/sqrt(ap_var)).^3,w');
    %ap_kurt = dotprod(((fapprox' - fhat(1,1))/sqrt(ap_var)).^4,w');
    fhat(1,1)
    error_mean = [error_mean;abs(fhat(1,1) - exact_mean)];
    %error_var = [error_var;abs(var(fapprox) - exact_variance)];
    %error_skew = [error_skew;abs(skewness(fapprox,1) - exact_skewness)];
    %error_kurt = [error_kurt;abs(kurtosis(fapprox,1) - exact_kurtosis)];
    %mu = [mu;fhat(1,1)];
    %var = [var;ap_var];
    %skew =[skew;ap_skew];
    %kurt = [kurt;ap_kurt];
    
end


%Plots for the 'UQ' stuff
figure(2)
mu = mean(u_sol,2);
sdev = std(u_sol')';
plot(x,mu+sdev,'b');
hold on;
plot(x,mu-sdev,'b');
x2 = [x, fliplr(x)];
inBetween = [mu+sdev, fliplr(mu-sdev)];
s = fill(x2, inBetween, 'b');
alpha(s, 0.5);
hold on;
plot(x,mu,'r');