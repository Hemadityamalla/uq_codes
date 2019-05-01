%This function has just one basis as the piecewise linear approximation of
%the integrand. Starting with a set of nodes (mother set), the accuracy of the approximation is kept fixed by taking
%an initial set of points and iteratively adding nodes from the mother set
%and implicitly generate a QR........
%run using the command (example): 
%[x,w,y] = cappedaccuracy_quad([5,5], @(x)exp(-x.^2/2)/sqrt(2*pi), 1e3);


function [x,w,y] = cappedaccuracy_quad(degree,f,Kmax)
D = degree(1); N = degree(2);
y = rand(Kmax,1); %Samples to be used
nodes = y(1:D); %The mother set

%Adding points and generating a quad rule
x = y(D+1:D+1+N-1);
w = ones(N,1)/(N);
iter = D+1+N;
node_iter = 1;
while (iter < Kmax-1 && node_iter <= D)
    %Adding a node from the mother set
    if node_iter <= D %this takes care of the finite size of the mother set
        x(end+1) = nodes(node_iter);
        w = [((iter-D-1+1)/(iter-D-1+2))*w;1./(iter-D-1+2)];
        node_iter = node_iter+1;
    end
    %Adding a node from the sample set
    x(end+1) = y(iter+1);
    %Rescaling the weights
    w = [((iter-D-1)/(iter+1-D-1))*w;1./(iter-D-1+1)];
    
    
    
    %Construct Vandermonde matrix
    %V = general_vandermonde(x, @(x,k) x.^(k-1), 1:N); %there is a possibility that I can go beyond N
    V = fliplr(vander(x))';
    
    %if node_iter>=3 %this makes sure that the function is not added until at least two nodes from the mother set are added
    %V(2,:) = interp1(nodes(1:node_iter-1),f(nodes(1:node_iter-1)),x,'linear','extrap');
    
    
    V(2,:) = interp1(nodes,f(nodes),x,'linear','extrap');
    
    
    %plot(0:0.01:1, interp1(nodes(1:node_iter-1),f(nodes(1:node_iter-1)),0:0.01:1,'linear','extrap'),'b-',sort(nodes(1:node_iter-1)),f(sort(nodes(1:node_iter-1))),'ro');
    %pause(0.1);
    %fprintf("Interpolation error: %f \n", norm(interp1(nodes(1:node_iter-1),f(nodes(1:node_iter-1)),0:0.01:1,'linear','extrap') - f(0:0.01:1)));
    %end
    
    %Applying Caratheodory
    nullVec = null(V(1:end-1,:)); %Neglecting the last row since we need a rectangular matrix(one basis less)
    c = nullVec(:,1);
    [alpha1, k1] = min(w(c>0)./c(c>0));
    [alpha2, k2] = max(w(c<0)./c(c<0));
    id1 = find(c > 0);
    id2 = find(c < 0);
    k1 = id1(k1);
    k2 = id2(k2);
    %check to determine which node to eliminate
    [~,newidx] = intersect(x,nodes,'stable');
    k = [k1;k2]; alpha = [alpha1;alpha2];
    comp_mat = (newidx' == k);
    ndel = sum(comp_mat,2);
    if isequal(ndel,[1;1])
        %there is no node removal here, moving on with an extra added node
        fprintf('iter: %d, no deletion \n',iter);
        iter = iter+1;
        continue;       
    elseif isequal(ndel,[1;0]) %k1 cannot be used
        k = k2; alpha = alpha2;
    elseif isequal(ndel,[0;1]) %k2 cannot be used
        k = k1; alpha = alpha1;
    else                        %Either k1,k2 can be used
        k = k1; alpha = alpha1;
    end
    w = w-alpha*c;
    x(k, :) = []; w(k) = [];
    iter = iter+1;
end
%     figure(1)
%     plot(0:0.01:1, interp1(nodes, f(nodes), 0:0.01:1,'linear'),'b.');
%     hold on;
%     plot(nodes, f(nodes),'ro');
%     figure(2);
%     xlim([0,1]);
%     title(num2str(N));
%     ylim([0,sum(degree)+1]);
%     scatter(x, N*ones(length(x),1),10,w,'filled');
%     hold on;
end
                                                                                                                                                                            
