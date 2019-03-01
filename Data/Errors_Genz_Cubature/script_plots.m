clear;clc;set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',16);
%Reading MC data
d = 6;

%Cinf fns
for ii=1:d
   Cinf.smolyak{ii} = csvread(strcat('./Cinf_fn/Smolyak_d',num2str(ii),'.dat'));
   Cinf.mc{ii} = reshape(csvread(strcat('./Cinf_fn/MC_d',num2str(ii),'.dat')),[d+1,2]);
   Cinf.tensor{ii} = csvread(strcat('./Cinf_fn/Tensor_d',num2str(ii),'.dat'));
   loglog(Cinf.smolyak{ii}(:,1), Cinf.smolyak{ii}(:,2), 'r');
   hold on;
   loglog(Cinf.mc{ii}(:,1), Cinf.mc{ii}(:,2),'b');
   hold on;
   loglog(Cinf.tensor{ii}(:,1), Cinf.tensor{ii}(:,2),'g');
end


%C0 fns
for ii=1:d
   C0.smolyak{ii} = csvread(strcat('./Cont_fn/Smolyak_d',num2str(ii),'.dat'));
   C0.mc{ii} = reshape(csvread(strcat('./Cont_fn/MC_d',num2str(ii),'.dat')),[d+1,2]);
   C0.tensor{ii} = csvread(strcat('./Cont_fn/Tensor_d',num2str(ii),'.dat'));
end

%Disc fns
for ii=1:d
   Disc.smolyak{ii} = csvread(strcat('./Disc_fn/Smolyak_d',num2str(ii),'.dat'));
   Disc.mc{ii} = reshape(csvread(strcat('./Disc_fn/MC_d',num2str(ii),'.dat')),[d+1,2]);
   Disc.tensor{ii} = csvread(strcat('./Disc_fn/Tensor_d',num2str(ii),'.dat'));
end
