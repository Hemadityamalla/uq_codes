clear;clc;close all;
x = 0.5*(1 + cos((0:100)*pi/100));
for n=1:10
   plot(x,haar(n,x));
   hold on;
end