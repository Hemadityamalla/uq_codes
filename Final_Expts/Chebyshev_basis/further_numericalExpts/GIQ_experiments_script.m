clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

rng(1,'twister'); %Seeding for reproducibility.


xpos = 500;ypos = 500; width = 1000; height = 800;

%Run the script for all the test cases and save the errors

for i=1:6 %functions
   for j=4 %quadrature
       fprintf("Function no. %i and Quadrature %i \n", i,j);
       Generalized_Implicit_Quad(j,i);
   end
end


%Load the scripts and plot the errors
N = 2:2:30;
for ifn = 1:6
    figure(ifn)
    for ires = 1:4
        data = dlmread(strcat('Errors_quad',num2str(ires),'_fn',num2str(ifn),'.dat'));
        loglog(N,data(:,1),'-o'); xlim([2,30]); 
        hold on;        
    end
    %xlabel('N');ylabel('Absolute error');legend('logarithmic','sinusoidal','polynomial');
    hold on;grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end