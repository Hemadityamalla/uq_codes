clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

%f = @(x) exp(-0.5*sum(abs(x-0.5),2));
f = @(x) cos(2*pi*0.3 + 0.5*(sum(x,2)));
%f = @(x) (prod(x,2) <= 0)*0 + (prod(x,2) > 0).*exp(0.5*sum(x,2)); %Discontinuous
d = 6; %dimension of the random vector

exact = (2^d)*cos(2*pi*0.3 + 0.25*d)*(sin(0.25)/0.5)^d;%0.5*(((1.0 - exp(-0.25))/0.5)^d + ((1.0 - exp(-0.25))/0.5)^d);%cos(2*pi*0.3 + 0.25*d)*(sin(0.25)/0.5)^d;%(2*(exp(0.5) - 1))^d;
error = [];
numpts = [];

for i=2*(10.^[1,2,3,4,5,6,7])
   eval_pts = rand(i,d);
   mu = mean(f(eval_pts))
   numpts(end+1) = numel(eval_pts);
   error(end+1) = abs(mu - exact);
end
csvwrite(strcat('MC_d',num2str(d),'.dat'),[numpts,error]);
polyfit(log(numpts),log(error),1)