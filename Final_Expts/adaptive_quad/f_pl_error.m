%Script to demonstrate the error in using piecewise linear approximation
%with linear extrapolation

clear;clc;close all;format long;


fn = genz_fns(0:0.001:1,rand(1,1),rand(1,1), 4);
Y = randn(1e6,1);
F = Y(1:20,1);
fpl = griddedInterpolant(sort(F),fn(sort(F)),'linear');
plot(sort(Y),fn(sort(Y)),'b',sort(F),fpl(sort(F)),'ro',sort(Y),fpl(sort(Y)),'y.')