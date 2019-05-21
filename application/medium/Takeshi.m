%clear;clc;format long;

%New wake model

%Velocity deficit
Ct = 0.36;
Ia = 3.5e-2;


k = @(Ct, Ia) 0.11*Ct.^1.07.*Ia.^0.2;
eps = @(Ct,Ia) 0.23*Ct.^(-0.25).*Ia.^0.17;




a = @(Ct,Ia) 1.93*Ct.^(-0.75).*Ia^0.17;
b = @(Ct,Ia) 0.42*Ct.^0.6.*Ia.^0.2;
c = @(Ct,Ia) 0.15*Ct.^(-0.25).*Ia.^(-0.7);


F = @(x,Ct,Ia) 1./(a(Ct,Ia) + b(Ct,Ia).*x + c(Ct,Ia).*(1 + x).^(-2)).^2;

x = 0:1:10;
r = 0;
sigma = k(Ct,Ia)*x + eps(Ct,Ia);
phi = exp(-r.^2./(2*sigma.^2));
du_uh = F(x,Ct,Ia).*phi;

plot(x,du_uh);xlim([0,10]);ylim([0,1]);
hold on;