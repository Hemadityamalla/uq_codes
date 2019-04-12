close all; clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',20);

qrule = {'clencurt','legendre'};
xpos = 500;ypos = 500; width = 1000; height = 800;

for testFn = 1:6
    figure(testFn);
for ii=1:length(qrule)
   fname = strcat('Error_',num2str(qrule{ii}),'_fn',num2str(testFn),'.dat');
   data = csvread(fname);
   semilogy(5*(1:length(data)),data,'-o');
   hold on;
end
for ii=1:6
   fname = strcat('Errors_quad',num2str(ii),'_fn',num2str(testFn),'.dat');
   data = csvread(fname);
   semilogy(5*(1:length(data)),data,'-*');
   hold on;
end

xlabel('Number of Points'); ylabel('Absolute Error'); xlim([4,44]);
legend('clencurt','legendre','log','p13','exp','sin','rat.','m13');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end