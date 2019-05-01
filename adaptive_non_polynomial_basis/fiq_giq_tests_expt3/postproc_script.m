clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

fns = {'gaussian','c0','disc'};
res = {'_005','_010','_025'};


%need to figure out how to customise a line type and a symbol for each plot
for ifn = 1:3
    figure(ifn)
    for ires = 1:3
        data = dlmread(strcat(fns{ifn},res{ires},'.dat'));
        
        semilogy(data(:,1),data(:,2),'r',data(:,1),data(:,3),'b');
        hold on;        
    end
    hold on;
end

figure(4)
for i=1:3
    data = dlmread(strcat(fns{i},res{1},'.dat'));
    semilogy(data(:,1),data(:,2),'r',data(:,1),data(:,3),'b');
    hold on;
end
