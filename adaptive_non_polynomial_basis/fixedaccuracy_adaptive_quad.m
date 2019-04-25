%This function has just one basis as the piecewise linear approximation of
%the integrand. The accuracy of the approximation is kept fixed by taking
%an initial set of points and iteratively adding nodes to this set to
%increase the accuracy of the QR upto Kmax.
%run using the command (example): 
%[x,w,y] = fixedaccuracy_adaptive_quad([5,5], @(x)exp(-x.^2/2)/sqrt(2*pi), 1e3)


function [x,w,y] = fixedaccuracy_adaptive_quad(degree,f,Kmax)
D = degree(1); N = degree(2);
y = rand(Kmax,1); %Samples to be used
nodes = y(1:D);%rand(D,1); %The fixed set of points

%Adding points and generating a quad rule
x = [nodes;y(D+1:N+D+1,1)];
w = ones(D+N,1)/(D+N);
iter = D+N+2;
replace = 0;
while iter < Kmax-1
    if replace == 0
        %Add a node and rescale the weights
        x = [x;y(iter,1)];
        w = ones(length(x),1)/length(x);%[(D+N)*w/(D+N+1);1./(D+N+1)];
    else
        newNode = rand(1,1);
        %y(iter) = newNode;
        x(end) = newNode; %Node replacement
        w = ones(length(x),1)/length(x);
    end
    
    %Construct Vandermonde matrix
    V = general_vandermonde(x, @(x,k) x.^(k-1), 1:length(nodes));
%     if length(nodes) >= 2
%        V(2,:) = interp1(nodes,f(nodes),x,'linear','extrap');
%     end
    
    %Applying Caratheodory
    nullVec = null(V);
    c = nullVec(:,1);
    [alpha1, k1] = min(w(c>0)./c(c>0));
    [alpha2, k2] = max(w(c<0)./c(c<0));
    id1 = find(c > 0);
    id2 = find(c < 0);
    k1 = id1(k1);
    k2 = id2(k2);
    %check to determine which node to eliminate
    [~,newidx] = setdiff(x,nodes);
    k = [k1;k2];alpha = [alpha1;alpha2];
    comp_mat = (newidx' == k);
    if sum(sum(comp_mat)) == 0
        %there is no node removal here, so replace the added sample
        %with a new sample
        replace = 1;
        fprintf('iter: %d, replacement iteration \n',iter);
        continue;
    else
        [kidx,~,~] = find(comp_mat,1,'last'); %gives the last occurence
        k = k(kidx); alpha = alpha(kidx);
        replace = 0;
        iter = iter + 1;
     end
    w = w-alpha*c;
    x(k, :) = []; w(k) = [];
end
    figure(1)
    plot(0:0.01:1, interp1(nodes, f(nodes), 0:0.01:1,'linear'),'b.');
    hold on;
    plot(nodes, f(nodes),'ro');
    figure(2);
    xlim([0,1]);
    title(num2str(N));
    ylim([0,sum(degree)+1]);
    scatter(x, N*ones(length(x),1),10,w,'filled');
    hold on;
end

