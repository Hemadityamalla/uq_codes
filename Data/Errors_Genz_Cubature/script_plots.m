clear;clc;set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',16);
%Reading MC data
d = 6;

%Cinf fns
for ii=1:d
   Cinf.smolyak{ii} = csvread(strcat('./Cinf_fn/Smolyak_d',num2str(ii),'.dat'));
   Cinf.mc{ii} = reshape(csvread(strcat('./Cinf_fn/MC_d',num2str(ii),'.dat')),[d+1,2]);
   Cinf.tensor{ii} = csvread(strcat('./Cinf_fn/Tensor_d',num2str(ii),'.dat'));
   %polyfit(log(Cinf.tensor{ii}(:,1)),log(Cinf.tensor{ii}(:,2)),1)
   %polyfit(log(Cinf.mc{ii}(:,1)),log(Cinf.mc{ii}(:,2)),1)
   %polyfit(log(Cinf.smolyak{ii}(:,1)),log(Cinf.smolyak{ii}(:,2)),1)
end


%C0 fns
for ii=1:d
   C0.smolyak{ii} = csvread(strcat('./Cont_fn/Smolyak_d',num2str(ii),'.dat'));
   C0.mc{ii} = reshape(csvread(strcat('./Cont_fn/MC_d',num2str(ii),'.dat')),[d+1,2]);
   C0.tensor{ii} = csvread(strcat('./Cont_fn/Tensor_d',num2str(ii),'.dat'));
   %polyfit(log(C0.tensor{ii}(:,1)),log(C0.tensor{ii}(:,2)),1)
   %polyfit(log(C0.mc{ii}(:,1)),log(C0.mc{ii}(:,2)),1)
   %polyfit(log(C0.smolyak{ii}(:,1)),log(C0.smolyak{ii}(:,2)),1)
end

%Disc fns
for ii=1:d
   Disc.smolyak{ii} = csvread(strcat('./Disc_fn/Smolyak_d',num2str(ii),'.dat'));
   Disc.mc{ii} = reshape(csvread(strcat('./Disc_fn/MC_d',num2str(ii),'.dat')),[d+1,2]);
   Disc.tensor{ii} = csvread(strcat('./Disc_fn/Tensor_d',num2str(ii),'.dat'));
   polyfit(log(Disc.tensor{ii}(:,1)),log(Disc.tensor{ii}(:,2)),1)
   %polyfit(log(Disc.mc{ii}(:,1)),log(Disc.mc{ii}(:,2)),1)
   %polyfit(log(Disc.smolyak{ii}(:,1)),log(Disc.smolyak{ii}(:,2)),1)
end




% %Plot for comparing the convergence of different methods at d=5
% xpos = 500;ypos = 500; width = 1000; height = 800;
% figure(1)
% loglog(Cinf.smolyak{5}(:,1),Cinf.smolyak{ii}(:,2),'bo-','MarkerFaceColor','b');
% hold on;
% loglog(Cinf.mc{5}(:,1),Cinf.mc{5}(:,2),'b s:','MarkerFaceColor','b');
% hold on;
% loglog(Cinf.tensor{5}(2:end,1),Cinf.tensor{5}(2:end,2),'b ^--','MarkerFaceColor','b');
% xlabel('Number of Points'); ylabel('Absolute Error'); 
% legend('Smolyak','Monte-Carlo','Tensor');
% grid on;set(gcf,'Position',[xpos ypos width height]); box on;
% figure(2)
% loglog(C0.smolyak{5}(:,1),C0.smolyak{ii}(:,2),'go-','MarkerFaceColor','g');
% hold on;
% loglog(C0.mc{5}(:,1),C0.mc{5}(:,2),'g s:','MarkerFaceColor','g');
% hold on;
% loglog(C0.tensor{5}(2:end,1),C0.tensor{5}(2:end,2),'g ^--','MarkerFaceColor','g');
% xlabel('Number of Points'); ylabel('Absolute Error');
% legend('Smolyak','Monte-Carlo','Tensor');
% grid on;set(gcf,'Position',[xpos ypos width height]); box on;
% figure(3)
% loglog(Disc.smolyak{5}(:,1),Disc.smolyak{ii}(:,2),'ro-','MarkerFaceColor','r');
% hold on;
% loglog(Disc.mc{5}(:,1),Disc.mc{5}(:,2),'r s:','MarkerFaceColor','r');
% hold on;
% loglog(Disc.tensor{5}(2:end,1),Disc.tensor{5}(2:end,2),'r ^--','MarkerFaceColor','r');
% grid on;set(gcf,'Position',[xpos ypos width height]); box on;
% xlabel('Number of Points'); ylabel('Absolute Error');
% legend('Smolyak','Monte-Carlo','Tensor');
% figure(4)
% loglog(Cinf.smolyak{5}(:,1),Cinf.smolyak{ii}(:,2),'bo-','MarkerFaceColor','b');
% hold on;
% loglog(C0.smolyak{5}(:,1),C0.smolyak{ii}(:,2),'go-','MarkerFaceColor','g');
% hold on;
% loglog(Disc.smolyak{5}(:,1),Disc.smolyak{ii}(:,2),'ro-','MarkerFaceColor','r');
% grid on;set(gcf,'Position',[xpos ypos width height]); box on;
% xlabel('Number of Points'); ylabel('Absolute Error');
% legend('C_\infty','C_0','Piecewise C_0');