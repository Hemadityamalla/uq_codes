%run using the command (example): 
%[x,w,y] = adaptive_quad(20, @(x)exp(-x.^2/2)/sqrt(2*pi), 1e3)

function [x,w,y] = adaptive_quad(degree,f,Kmax)
nodes = 0.5;
increments = 2;
for N=2:increments:degree
    %at any instant length of nodes is (N-1) (if incrementing N by 1)
    y = rand(Kmax,1);
    %Initialize quad rule
    x = [nodes; y(1:(N-length(nodes)))]; %essentially you are adding only one node here(if incrementing N by 1)
    %x = y(1:N);
    w = ones(N,1)/(N);
    %Implicit quad rule
    replace = 0; D = N;
    
    while (D <= (Kmax-1))
       
       if replace == 0
        x = [x;y(D+1)]; %Node addition
        w = [((D)/(D+1))*w;1./(D+1)];
       else
           newNode = rand(1,1);
           y(D+1) = newNode;
          x(end) = newNode; %Node replacement
          w = ones(length(x),1)/length(x);
       end
       
       
       V = zeros(N, length(x));
       V(1,:) = ones(1,length(x));
       i=1;
       for ii=2:length(nodes)
           V(ii,:) = interp1(nodes(1:ii),f(nodes(1:ii)),x,'linear','extrap'); 
           i=i+1;
       end
       %Filling up the rest of the rows with polynomials+
%        for jj=i+1:N %Must'nt this be i:D
%           V(jj,:) = x.^(jj-1);
%        end
%--------------Expt by adding low degree polynomials, OG code above
        iter = 2;
       for jj=i+1:N
          V(jj,:) = x.^(iter-1);
          iter = iter+1;
       end
       nullVec = null(V);
       c = nullVec(:,1);

       % Apply Caratheodorys' and determine the node that will be removed
        [alpha1, k1] = min(w(c>0)./c(c>0));
        [alpha2, k2] = max(w(c<0)./c(c<0));
        id1 = find(c > 0);
        id2 = find(c < 0);
        k1 = id1(k1);
        k2 = id2(k2);
        %check to determine which node to eliminate
        newidx = length(x(1:end-increments)):length(x);
        k = [k1;k2];alpha = [alpha1;alpha2];
        comp_mat = (newidx == k);
        if sum(sum(comp_mat)) == 0
            %there is no node removal here, so replace the added sample
            %with a new sample
            replace = 1;
            %fprintf('N: %d, D = %d, replacement iteration \n',N,D);
            continue;
        else
            [kidx,~,~] = find(comp_mat,1,'last'); %gives the last occurence
            k = k(kidx); alpha = alpha(kidx);
            replace = 0;
            D = D + 1;
         end


%     %Update weights
%      [alpha, k] = min(w(c>0)./c(c>0));
% 
%         id = find(c > 0);
% 
%         k = id(k);
%  
    w = w-alpha*c;
    x(k, :) = []; w(k) = [];
%         
        
    end
    %fprintf("Length of x: %d, Length of nodes: %d \n", length(x), length(nodes));
    nodes = x;
%     figure(1)
%     plot(0:0.01:1, interp1(nodes, f(nodes), 0:0.01:1,'linear'),'b.');
%     hold on;
%     plot(nodes, f(nodes),'ro');
%     figure(2);
%     xlim([0,1]);
%     title(num2str(N));
%     ylim([0,degree+1]);
%     scatter(nodes, N*ones(length(nodes),1),10,w,'filled');
%     hold on;
end

end

