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
    for ires = [1,2,3,5,6]
        data = dlmread(strcat('Errors_quad',num2str(ires),'_fn',num2str(ifn),'.dat'));
        if ires == 5 || ires == 6
            loglog(data(:,1), data(:,2),'-o'); xlim([2,30]);
        else
            loglog(N,data(:,1),'-o'); xlim([2,30]);
        end
        hold on;        
    end
    xlabel('N');ylabel('Absolute error');legend('\beta_1','\beta_2','monomial (\beta_3)','\beta_4','\beta_5');
    hold on;grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end