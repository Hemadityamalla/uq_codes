clear;clc; format long;


%points used for piecewise linear interpolation
Kmax = 1e5;

f = @(x) 1.0*(x > 0.5);%Test function to be integrated

exact = 0.5;
error_adaptive = [];
error_trad = [];
D = 14;
for degree = 4:2:14
    %[xa,wa] = adaptive_quad(degree,f,Kmax);
    [xt,wt] = fixed_implict_quad(degree,Kmax);
    %error_adaptive = [error_adaptive;abs(dotprod(f(xa'),wa)-exact)];
    %error_trad = [error_trad;abs(exact - dotprod(f(xt'),wt))];
    figure(2);
    %xlim([0,1]);

    ylim([0,D+1]);
    scatter(xt, degree*ones(length(xt),1),10,wt,'filled');
    hold on;
end
%figure(2);
%semilogy(5:5:60,error_adaptive,'b',5:5:60,error_trad,'r');
