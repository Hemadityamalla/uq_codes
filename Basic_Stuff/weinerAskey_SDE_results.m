%clear;clc;

t0 = 0.0; tf = 1.0;

ic = 1.0;

k = randn(1,1); 
[t,y] = ode45(@(t,y) -k*y, [t0,tf],ic);
plot(t,y);

