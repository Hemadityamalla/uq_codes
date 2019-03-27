close all; clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',20);

qrule = {'clencurt','log','m13','p13','legendre'};
xpos = 500;ypos = 500; width = 1000; height = 800;

for testFn = 1:6
    figure(testFn);
for ii=1:length(qrule)
   fname = strcat('Error_',num2str(qrule{ii}),'_fn',num2str(testFn),'.dat');
   data = csvread(fname);
   semilogy(5*(1:length(data)),data,'-o');
   hold on;
end

xlabel('Number of Points'); ylabel('Absolute Error'); xlim([4,44]);
legend(qrule);
grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end