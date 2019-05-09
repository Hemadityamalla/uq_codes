%This function has just one basis as the piecewise linear approximation of
%the integrand. Starting with a set  fixedNodes (mother set), the accuracy of the approximation is kept fixed by taking
%an initial set of points and iteratively adding nodes from the mother set
%and implicitly generate a QR........
%run using the commands (example):
%N = 10; D = 5;[x,w,fixedNodes,y] = cappedaccuracy_quad_v2([D,N], @(x)exp(-x.^2/2)/sqrt(2*pi), 5e2);


%D -- fixed number of nodes that are used to construct a piecewise linear
%approximation f^D(x)
%N -- Exactly integrates polynomials of degree upto (N-2) and f^D(x). N is
%preferrably at least 2*D.

%The following commands can be ran to validate the QR:
% N = 20;D = 5;[x,w,fixedNodes,y] = cappedaccuracy_quad_v2([D,N], @(x)exp(-x.^2/2)/sqrt(2*pi), 5e2);
%sum(w) %Must be close to 1.0
%sum(x.*w) == mean(y)
%sum(x.^(N-2).*w) == mean(y.^(N-2))
%pp = griddedInterpolant(sort(fixedNodes), f(sort(fixedNodes)),'linear');
%sum(pp(x).*w) == mean(pp(x))


function [x,w,nodes,Y] = cappedaccuracy_quad_v2(degree,f,Kmax)
D = degree(1); N = degree(2);
Y = rand(Kmax,1); %Samples to be used
nodes = Y(1:D); %Fixed sample set
y = Y(D+1:end); leftOverSize = length(y);

%Initialize starting QR
x = y(1:N);
w = ones(N,1)/(N);
iter = N;
leftOverIter = N;
added_nodes = 0;
%Constructing the piecewise linear interpolant
pp = griddedInterpolant(sort(nodes), f(sort(nodes)),'linear');
while (iter < leftOverSize+D)
%         figure(1)
%     %scatter(nodes, (iter-1)*ones(length(nodes),1),10,'filled');
%     hold on;
%     scatter(x, iter*ones(length(x),1),10,w,'filled');
    if added_nodes < D
        x(end+1) = nodes(added_nodes+1); %Adding a node from the fixed sample set
        added_nodes = added_nodes+1;
    elseif added_nodes >= D
        x(end+1) = y(leftOverIter+1);
        leftOverIter = leftOverIter + 1;
    end
    %Rescaling the weights
    %initL = length(w);
    w = [(iter/(iter+1))*w;1/(iter+1)];

    %Construct Vandermonde matrix
    V = general_vandermonde(x, @(x,k) x.^(k-1), 1:N); 
    
    
    V(end,:) = pp(x);
    
    
    %plot(0:0.01:1, interp1(nodes(1:node_iter-1),f(nodes(1:node_iter-1)),0:0.01:1,'linear','extrap'),'b-',sort(nodes(1:node_iter-1)),f(sort(nodes(1:node_iter-1))),'ro');
    %pause(0.1);
    %fprintf("Interpolation error: %f \n", norm(interp1(nodes(1:node_iter-1),f(nodes(1:node_iter-1)),0:0.01:1,'linear','extrap') - f(0:0.01:1)));
    %end
    
    %Applying Caratheodory
    nullVec = null(V);
    c = nullVec(:,1);
    if size(nullVec,2) > 1
        size(nullVec,2)
    end
    [alpha1, k1] = min(w(c>0)./c(c>0));
    [alpha2, k2] = max(w(c<0)./c(c<0));
    id1 = find(c > 0);
    id2 = find(c < 0);
    k1 = id1(k1);
    k2 = id2(k2);
    %check to determine which node to eliminate
    [~,newidx] = intersect(x,nodes,'stable'); %Gives only the last element
    %newidx
    k = [k1;k2]; alpha = [alpha1;alpha2];
    comp_mat = (newidx' == k);
    ndel = sum(comp_mat,2);
    if isequal(ndel,[1;1]) %This step here is causing more nodes to be added and nothing to be removed
        %there is no node removal here, moving on with an extra added node
        %fprintf('iter: %d, no deletion \n',iter); 
        iter = iter+1;
        continue;       
    elseif isequal(ndel,[1;0]) %k1 cannot be used
        k = k2; alpha = alpha2;
    elseif isequal(ndel,[0;1]) %k2 cannot be used
        k = k1; alpha = alpha1;
    else                        %Either k1,k2 can be used
        %fprintf("Two nodes can be deleted\n");
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
                                                                                                                                                                            
