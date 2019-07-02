clear;clc; close all; format long
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

N = 5; D = 2; level = 5; Kmax = 1e5;

for fnType = 1:6
    samples = rand(Kmax,1);
   figure(fnType)
   [QR] = cappedaccuracy_quad_nested([D,N], level,genz_fns(0:0.001:1, rand(1,1), rand(1,1), fnType), samples);
end