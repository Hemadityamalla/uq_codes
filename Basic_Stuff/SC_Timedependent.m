clear;
clc;

N = 50; l = 1.0;
x = linspace(0,l,N);


dx = l/N;
dt = 0.00001;
tfinal = 0.1;
To = 300;



for i=[1:1:10e3]
    Tinit = To + To*sin(pi*x(2:N-1))';
    T = zeros(size(Tinit));
    b = zeros(N-2,1);
    alpha = exp(1.0 + 0.75*randn); %Log-normal distn- always gives a positive sample
    K = alpha*dt/(dx^2);
    b(1,1) = K*To; b(end,1) = K*To;
    A = gallery('tridiag',N-2,K, 1 - 2*K, K);
    t = 0;
    while t < tfinal
        T = A*Tinit + b;
        t = t + dt;
        %plot(x,[To;T;To],'r-.');
        %hold on;
        %pause(0.001);
        Tinit = T;
    end
    Tsample(i,1) = T(17,1);
end
Tpos = Tsample(Tsample > 0);
histogram(Tpos(Tpos < 1000));