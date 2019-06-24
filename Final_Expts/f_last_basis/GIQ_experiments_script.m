clear;clc;close all;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

%Run the script for all the test cases and save the errors

for i=1:6 %functions
   for j=1:3 %quadrature
       fprintf("Function no. %i and Quadrature %i \n", i,j);
       Generalized_Implicit_Quad(j,i);
   end
end


%Load the scripts and plot the errors
N = 2:2:30;
for ifn = 1:6
    figure(ifn)
    for ires = 1:3
        data = dlmread(strcat('Errors_quad',num2str(ires),'_fn',num2str(ifn),'.dat'));
        semilogy(N,data(:,1),'-o'); xlim([2,30]); 
        hold on;        
    end
    legend('log','sin','polyn.');
    hold on;
end