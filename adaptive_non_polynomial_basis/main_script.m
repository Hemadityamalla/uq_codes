clear;clc; format long;


%points used for piecewise linear interpolation
Kmax = 1e4;

%f = @(x) 1.0*(x > 0.5);%Test function to be integrated
f = @(x) x.*(x < 0.5) + (x >= 0.5).*exp(-0.5*x);
%f = @(x) exp(-x.^2/2)/sqrt(2*pi);

exact = 0.5;
error_adaptive = [];
error_trad = [];
D = 55;
range = 5:10:D;
for degree = range
    [xa,wa,ya] = almost_adaptive_quad(degree,f,Kmax);
    [xt,wt,yt] = fixed_implict_quad(degree,Kmax);
    error_adaptive = [error_adaptive;abs(dotprod(f(xa'),wa)-mean(f(ya)))];
    error_trad = [error_trad;abs(mean(f(yt)) - dotprod(f(xt'),wt))];
    figure(1);
    %xlim([0,1]);

    ylim([0,D+1]);
    scatter(xt, degree*ones(length(xt),1),10,wt,'filled');
    hold on;
end
figure(2);
semilogy(range,error_adaptive,'b',range,error_trad,'r');
