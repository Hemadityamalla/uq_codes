clear;clc;close all;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);


%Quadratures
figure(1)
for i=1:20
   lim = 1;
   %[x,w] = quadgen(i,'ClenshawCurtis');
   % addpath('/home/hemaditya/Documents/Approximation_Theory/chebfun')
   [x,~]= h(i);%(-cos((pi*((1:i)-1))/(i-1)) + 1)*0.5;
   x = 0.5*(x+1);
  
   %if  (floor(log(i-1)/log(2)) ==  log(i-1)/log(2)) && (i ~= 2)
        %scatter(x,length(x)*ones(length(x),1),25,'ro','filled');
   %else
       scatter(x,length(x)*ones(length(x),1),25,'bo','filled');
   %end
   hold on;
   scatter(-lim:0.005:lim,i*ones(length(-lim:0.005:lim),1),5,'k.');
   set(gcf,'Position',[0 0 800 800]); box on;
   %cbh = colorbar('northoutside');%,'Ticks',[0 0.2 0.4 0.6 0.8 1],...
       %'TickLabels',{'0.0','0.2','0.4','0.6','0.8','1'});
       %set(cbh,'YTick',[0:0.2:2])
   xlim([0,lim]);ylim([0,20]);
   xlabel('\xi');ylabel('N')
   hold on;
end

%Cubatures
figure(2)
lim = 1;d=2;k=5;
[x,~] = smolyakSparseGrid(d,k,@(x) 2^(x-1), 'ClenshawCurtis');
x = 0.5*(x+1);
scatter(x(1,:),x(2,:),25,'bo','filled');
%hold on;
%scatter(-lim:0.005:lim,i*ones(length(-lim:0.005:lim),1),5,'k.');
set(gcf,'Position',[0 0 800 800]); box on;
%colorbar('northoutside','Ticks',[0 0.2 0.4 0.6 0.8 1],...
 %  'TickLabels',{'0.0','0.2','0.4','0.6','0.8','1'});
xlim([0,lim]);ylim([0,lim]);
xlabel('\xi_1');ylabel('\xi_2')
%hold on;
figure(5)
[xi,w] = quadgen(2^(k),'ClenshawCurtis');
x = setprod(xi,d);
x = 0.5*(x+1);
w = prod(setprod(w,d),2);
scatter(x(:,1),x(:,2),25,'bo','filled');
%hold on;
%scatter(-lim:0.005:lim,i*ones(length(-lim:0.005:lim),1),5,'k.');
set(gcf,'Position',[0 0 800 800]); box on;
%colorbar('northoutside','Ticks',[0 0.2 0.4 0.6 0.8 1],...
 %  'TickLabels',{'0.0','0.2','0.4','0.6','0.8','1'});
xlim([0,lim]);ylim([0,lim]);
xlabel('\xi_1');ylabel('\xi_2')
%hold on;


%Reduced quadrature
figure(3)
[x,w] = quadgen(20,'Hermite');
xmax = max(x);
%x = 0.5*(x+1);
for i=1:20
   lim = xmax+0.5;
   scatter(x,length(x)*ones(length(x),1),25,'bo','filled');
   hold on;
   scatter(-lim:(lim/100):lim,i*ones(length(-lim:(lim/100):lim),1),5,'k.');
   set(gcf,'Position',[0 0 800 800]); box on;
   %colorbar('northoutside','Ticks',[0 0.2 0.4 0.6 0.8 1],...
    %   'TickLabels',{'0.0','0.2','0.4','0.6','0.8','1'});
   xlim([-lim,lim]);
   ylim([0,20]);
   xlabel('\xi');ylabel('N')
   hold on;
   [x,w] = reduced_quad(x,w);
end

%Implicit_quad
figure(4)
Kmax = 5e3;
Y = randn(Kmax,1);
xmax = max(Y);
for i=2:20
   lim = xmax + 0.5;
   
   [x,w] = fixed_implict_quad(i,Y);
   %x = 0.5*(x+1);
   scatter(x,length(x)*ones(length(x),1),25,'bo','filled');
   hold on;
   scatter(-lim:(lim/100):lim,i*ones(length(-lim:(lim/100):lim),1),5,'k.');
   set(gcf,'Position',[0 0 800 800]); box on;
   %colorbar('northoutside','Ticks',[0 0.2 0.4 0.6 0.8 1],...
   %    'TickLabels',{'0.0','0.2','0.4','0.6','0.8','1'});
   xlim([-lim,lim]);ylim([0,20]);
   xlabel('\xi');ylabel('N')
   hold on;
end


