close all; clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',16);

f = @(x,t) 1./(1 + t*x);
n = 50;
t1 = linspace(-0.8,0.8,n); %Equispaced parameters
t2 = -1 + 2*rand(n,1); %Uniform random distribution (no guarantee that the t's are pairwise unique)
x = -1:0.01:1;
fn = zeros(n, length(x));
for ii=1:n
    fn(ii,:) = f(x,t1(ii));
    plot(x,fn(ii,:));
    hold on;
    pause(0.1);
end

[q,~,~] = qr(fn');
