close all; clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',20);

addpath('/ufs/hemadity/Documents/chebfun');

%Loading the quadrature rules

% %Logarithmic
% qrule = 'log';
% for ii=1:8
%    ggq.q{ii} = csvread(strcat('ggq_ln_',num2str(5*ii),'.csv'));
% end

% %x^(1/3)
% qrule = 'p13';
% for ii=1:4
%    ggq.q{ii} = (dlmread(strcat('ggq_1b3_',num2str(5*ii),'.dat')))';
% end
% 
% %x^(-1/3)
% qrule = 'm13';
% for ii=1:4
%    ggq.q{ii} = (dlmread(strcat('ggq_m13_',num2str(5*ii),'.dat')))';
% end

qrule = 'clencurt';
for ii=1:8
    [x,w] = clencurt(5*ii - 1);
    ggq.q{ii} = [0.5*(x+1),w*0.5];
end

% qrule = 'legendre';
% for ii=1:8
%    [x,w] = legpts(5*ii);
%    ggq.q{ii} = [0.5*(x+1),w'*0.5];
% end

for testFn = 1:6
    switch testFn
        case 1
            %Oscillatory
            u = @(x) cos(2*pi*0.6 + 0.5*x);
            exact = 4*cos(0.25 + 2*pi*0.6)*sin(0.25);
        case 2
            %product peak
            u = @(x) (0.5^(-2) + x.^2).^(-1);
            exact = 0.5*atan(0.5);
        case 3
            %Corner peak
            u = @(x) 1./(1 + 5*x).^2;
            exact = 1/6;
        case 4
            %C0
            u = @(x) 0.25*exp(-0.5*abs(x-0.5));
            exact = (1 - exp(-0.25));
        case 5
            %simple monomial
            u = @(x) x.^9;
            exact = 0.1;
        case 6
            %Discontinuous function
            u = @(x) (x <= 0.5)*0 + (x > 0.5).*exp(0.5*x);
            exact = 2*exp(0.5) - 2*exp(0.25);
    end
    fname = strcat('Error_',num2str(qrule),'_fn',num2str(testFn),'.dat');
    err = 0;
for ii=1:length(ggq.q)
    x = ggq.q{ii}(:,1);
    w = ggq.q{ii}(:,2);
    err = abs(exact - dotprod(u(x)',w));
    if ii==1
        dlmwrite(fname,err);
    else
        dlmwrite(fname,err,'-append');
    end
    
    
end
t = 10;
end




