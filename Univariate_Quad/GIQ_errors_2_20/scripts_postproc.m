close all;clear;clc;set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',16); format long;

nquad = 4;
N = 2:2:20;
%Loading the functions
iter = 1;
for ii = [1:nquad]
    %F1
    error.f1{iter} = [N',csvread(strcat('Errors_quad',num2str(ii),'_fn1.dat'))];
    %F2
    error.f2{iter} = [N',csvread(strcat('Errors_quad',num2str(ii),'_fn2.dat'))];
    %F3
    error.f3{iter} = [N',csvread(strcat('Errors_quad',num2str(ii),'_fn3.dat'))];
    %F4
    error.f4{iter} = [N',csvread(strcat('Errors_quad',num2str(ii),'_fn4.dat'))];
    %F5
    error.f5{iter} = [N',csvread(strcat('Errors_quad',num2str(ii),'_fn5.dat'))];
    iter = iter + 1;
end

% error_f1 = [];
% error_f2 = [];
% error_f3 = [];
% error_f4 = [];
% error_f5 = [];
% 
% %Loading ggq_ln
% for ii=1:8
%    ggq_ln.q{ii} = csvread(strcat('ggq_ln_',num2str(5*ii),'.csv'));
%    x = ggq_ln.q{ii}(:,1);
%    w = ggq_ln.q{ii}(:,2);
%     %1) Integration with monomial of order N-1
%     exact = 1/(5*ii);
%     approx=dotprod(x'.^(5*ii-1),w);
%     error_f1(end+1) = abs(exact - approx);
% 
%     %2) Integration with cont. Genz function
%     exact = 1 - exp(-0.25);
%     approx = dotprod(0.25*exp(-0.5*abs(x' - 0.5)),w);
%     error_f2(end+1) = abs(exact - approx);
%     
%     %4) Integration with corner peak Genz function
%     exact = 1./6;
%     approx = dotprod(1./(1 + 5*x').^2, w);
%     error_f3(end+1) = abs(exact - approx);
%     nQuads
%     %5) Integration with corner peak Genz function
%     exact = -1./2116;
%     approx = dotprod((x'.^45).*log(x'), w);
%     error_f4(end+1) = abs(exact - approx);
%     
%     %6) Integration with highly oscillatory function
%     exact = (2./55)*sin(55./2)^2;
%     approx = dotprod(sin(55*x'),w);
%     error_f5(end+1) = abs(exact - approx);
% end
% 
% error.f1{end+1} = [(5:5:40)',error_f1'];
% error.f2{end+1} = [(5:5:40)',error_f2'];
% error.f3{end+1} = [(5:5:40)',error_f3'];
% error.f4{end+1} = [(5:5:40)',error_f4'];
% error.f5{end+1} = [(5:5:40)',error_f5'];
% 
% error_f1 = [];nQuads
% error_f2 = [];
% error_f3 = [];
% error_f4 = [];
% error_f5 = [];
% 
% %Loading ggq_1b3
% for ii=1:4
%     ggq_1b3.q{ii} = (dlmread(strcat('ggq_1b3_',num2str(5*ii),'.dat')))';
%    x = ggq_ln.q{ii}(:,1);
%    w = ggq_ln.q{ii}(:,2);
%     %1) Integration with monomial of order N-1
%     exact = 1/(5*ii);
%     approx=dotprod(x'.^(5*ii-1),w);
%     error_f1(end+1) = abs(exact - approx);
% 
%     %2) Integration with cont. Genz function
%     exact = 1 - exp(-0.25);
%     approx = dotprod(0.25*exp(-0.5*abs(x' - 0.5)),w);
%     error_f2(end+1) = abs(exact - approx);
%     
%     %4) Integration with corner peak Genz function
%     exact = 1./6;
%     approx = dotprod(1./(1 + 5*x').^2, w);
%     error_f3(end+1) = abs(exact - approx);
%     
%     %5) Integration with corner peak Genz function
%     exact = -1./2116;
%     approx = dotprod((x'.^45).*log(x'), w);
%     error_f4(end+1) = abs(exact - approx);
%     
%     %6) Integration with highly oscillatory function
%     exact = (2./55)*sin(55./2)^2;
%     approx = dotprod(sin(55*x'),w);
%     error_f5(end+1) = abs(exact - approx);
% end
% 
% error.f1{end+1} = [(5:5:20)',error_f1'];
% error.f2{end+1} = [(5:5:20)',error_f2'];
% error.f3{end+1} = [(5:5:20)',error_f3'];
% error.f4{end+1} = [(5:5:20)',error_f4'];
% error.f5{end+1} = [(5:5:20)',error_f5'];

%Plots
xpos = 500;ypos = 500; width = 1000; height = 800;nQuads = 7; xrange = [1,25];
figure(1)
for ii=1:nquad
   semilogy(error.f1{ii}(:,1), error.f1{ii}(:,2),'-o'); 
   hold on;
end
xlabel('Number of Points'); ylabel('Absolute Error'); xlim(xrange);
legend('Q1','Q2','Q3','Q4','Q5','Q1 lit.','Q2 lit.');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;


figure(2)
for ii=1:nquad
   semilogy(error.f2{ii}(:,1), error.f2{ii}(:,2),'-o'); 
   hold on;
end
xlabel('Number of Points'); ylabel('Absolute Error'); xlim(xrange);
legend('Q1','Q2','Q3','Q4','Q5','Q1 lit.','Q2 lit.');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;

figure(3)
for ii=1:nquad
   semilogy(error.f3{ii}(:,1), error.f3{ii}(:,2),'-o'); 
   hold on;
end
xlabel('Number of Points'); ylabel('Absolute Error'); xlim(xrange);
legend('Q1','Q2','Q3','Q4','Q5','Q1 lit.','Q2 lit.');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;

figure(4)
for ii=1:nquad
   semilogy(error.f4{ii}(:,1), error.f4{ii}(:,2),'-o'); 
   hold on;
end
xlabel('Number of Points'); ylabel('Absolute Error'); xlim(xrange);
legend('Q1','Q2','Q3','Q4','Q5','Q1 lit.','Q2 lit.');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;


figure(5)
for ii=1:nquad
   semilogy(error.f5{ii}(:,1), error.f5{ii}(:,2),'-o'); 
   hold on;
end
xlabel('Number of Points'); ylabel('Absolute Error'); xlim(xrange);
legend('Q1','Q2','Q3','Q4','Q5','Q1 lit.','Q2 lit.');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;
