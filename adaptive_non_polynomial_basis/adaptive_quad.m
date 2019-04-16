
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
    for D=N:(Kmax-1)
       %Node addition
       x = [x;y(D+1)];
       w = [((D)/(D+1))*w;1./(D+1)];
       %Update weights
       V = zeros(N, length(x));
       V(1,:) = ones(1,length(x));
       i=1;
       for ii=2:length(nodes)
           V(ii,:) = interp1(nodes(1:ii),f(nodes(1:ii)),x,'linear','extrap'); 
           i=i+1;
       end
       %for jj=i:N %Must'nt this be i:D
       %    V(jj,:) = x.^(jj-1);  
       %end

       nullVec = null(V);
       c = nullVec(:,1);

       % Apply Caratheodory's and determine the node that will be removed
        [alpha, k] = min(w(c>0)./c(c>0));
        %[alpha, k] = max(w(c<0)./c(c<0));
        id = find(c > 0);%find(c < 0)
        k = id(k);

        %Node removal
        w = w-alpha*c;
        x(k, :) = []; w(k) = [];
        
    end
    %fprintf("Length of x: %d, Length of nodes: %d \n", length(x), length(nodes));
    nodes = x;
    plot(0:0.01:1, interp1(nodes, f(nodes), 0:0.01:1,'linear'),'.');
    %scatter(nodes, N*ones(length(nodes),1),10,w,'filled');
    %hold on;
    xlim([0,1]);
    pause(0.01);
    %ylim([0,degree+1])

end

end

