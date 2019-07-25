%Script solving the stochastic, steady state burgers' equation

clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);
addpath('/ufs/hemadity/Documents/chebfun');
xpos = 500;ypos = 500; width = 1000; height = 800;
rng(1,'twister');


xmesh = linspace(-1,1,20);
solinit = bvpinit(xmesh, @guess);

soln = bvp4c(@bvpfcn, @bcfcn, solinit);



plot(soln.x, soln.y(1,:), '-o');





function dydx = bvpfcn(x,y,nu)
    %nu = 0.05;
    dydx = zeros(2,1);
    dydx = [y(2), (y(1)*y(2))/nu];

end

function res = bcfcn(ya, yb)
    res = [ya(1)-1, yb(1)+1];
end

function g = guess(x)
    g = [sin(x),cos(x)];
end