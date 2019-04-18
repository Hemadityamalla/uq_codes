
%run using the command: 
%[x,w] = adaptive_quad(20, @(x)exp(-x.^2/2)/sqrt(2*pi), 1e3)

function [x,w] = adaptive_quad(degree,f,Kmax)
nodes = 0.5;
for N=2:degree
    y = rand(Kmax,1);
    %Initialize quad rule
    x = [nodes; y(1:(N-length(nodes)))];
    w = ones(N,1)/(N);
    %Implicit quad rule
    replace = 0;
    for D=N:(Kmax-1)
       
       if replace == 0
        x = [x;y(D+1)]; %Node addition
        w = [((D)/(D+1))*w;1./(D+1)];
       else
          x(end) = rand(1,1); %Node replacement
          w = ones(length(x),1)/length(x);
       end
       
       %Update weights
       V = zeros(N, length(x));
       V(1,:) = ones(1,length(x));
       i=1;
       for ii=2:length(nodes)
           V(ii,:) = interp1(nodes(1:ii),f(nodes(1:ii)),x,'linear','extrap'); 
           i=i+1;
       end
       for jj=ii+1:N %Must'nt this be i:D
          V(jj,:) = x.^(jj-1);  
       end

       nullVec = null(V);
       c = nullVec(:,1);

       % Apply Caratheodorys' and determine the node that will be removed
%         [alpha1, k1] = min(w(c>0)./c(c>0));
%         [alpha2, k2] = max(w(c<0)./c(c<0));
%         id1 = find(c > 0);
%         id2 = find(c < 0);
%         k1 = id1(k1);
%         k2 = id2(k2);
        %check to determine which node to eliminate
%         newidx = length(x(1:end-1)):length(x);
%         if sum(sum(newidx == [k1,k2]')) == 0
%             %there is no node removal here, so replace the added sample
%             %with a new sample
%             replace = 1;
%             continue;
%         else
%             k = k2; alpha = alpha2; %This is chosen arbitrarily--there is scope for better decision making
%             %Node removal (two possibilites)
%             
%             replace = 0;
%         end

     [alpha, k] = min(w(c>0)./c(c>0));

        id = find(c > 0);

        k = id(k);
 
    w = w-alpha*c;
    x(k, :) = []; w(k) = [];
        
        
    end
    %fprintf("Length of x: %d, Length of nodes: %d \n", length(x), length(nodes));
    nodes = x;
    figure(1)
    plot(0:0.01:1, interp1(nodes, f(nodes), 0:0.01:1,'linear'),'b.');
    hold on;
    plot(nodes, f(nodes),'ro');
    figure(2);
    xlim([0,1]);
    title(num2str(N));
    ylim([0,degree+1]);
    scatter(nodes, N*ones(length(nodes),1),10,w,'filled');
    hold on;
end

end

