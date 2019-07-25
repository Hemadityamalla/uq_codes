%Script to compare fpl with the best approximation polynomial of a certain
%degree

%clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');



fn = genz_fns(0:0.001:1,rand(1,1),rand(1,1), 5);
Y = rand(1e6,1); D = 10;
F = Y(1:D);

%Constructing piecewise linear approximant
fpl = griddedInterpolant(sort(F),fn(sort(F)'),'linear');

%Constructing the best approximation of degree N-2
error = [];
for N = 3:20
fbest = minimax(fn, N-2);

error(end+1) = norm(fbest(sort(Y)) - fn(sort(Y)),inf);
end
%plot(sort(Y), fbest(sort(Y)), 'bo', sort(Y), fpl(sort(Y)),'r-', sort(Y), fn(sort(Y)),'g.');
loglog(3:20, error,'bo-')
