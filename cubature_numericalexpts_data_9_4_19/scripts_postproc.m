close all; clear;clc;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',15);

qrule = {'MC_','Smolyak_','Smolyak_ggq_','Smolyak_RClenshaw_','Tensor_'};
xpos = 500;ypos = 500; width = 1000; height = 800;

for testFn = 1:6
    figure(testFn);
for ii=1:length(qrule)
   fname = strcat(qrule{ii},'f',num2str(testFn),'d5','.dat');
   data = dlmread(fname);
   loglog(data(:,1),data(:,2),'-o');
   hold on;
end
% for ii=1:6
%    fname = strcat('Errors_quad',num2str(ii),'_fn',num2str(testFn),'.dat');
%    data = csvread(fname);
%    semilogy(5*(1:length(data)),data,'-*');
%    hold on;
% end

xlabel('Number of Points'); ylabel('Absolute Error');xlim([1,1e6]);
legend('Monte-Carlo','Smolyak+CC','Smolyak+GGQ','Smolyak+CC.red.','Tensor');
grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end