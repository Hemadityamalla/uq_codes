%This function has just one basis as the piecewise linear approximation of
%the integrand. Starting with a set  fixedNodes (mother set), the accuracy of the approximation is kept fixed by taking
%an initial set of points and iteratively adding nodes from the mother set
%and implicitly generate a QR........
%run using the commands (example):
%N = 10; D = 5;[x,w,fixedNodes,y] = cappedaccuracy_quad_v3([D,N], @(x)exp(-x.^2/2)/sqrt(2*pi), 5e2);


%D -- fixed number of nodes that are used to construct a piecewise linear
%approximation f^D(x)
%N -- Exactly integrates polynomials of degree upto (N-2) and f^D(x). N is
%preferrably at least 2*D.

%The following commands can be ran to validate the QR:
% N = 20;D = 5;y=rand(5e2,1);[x,w,fixednodes] = cappedaccuracy_quad_v3([D,N], @(x)exp(-x.^2/2)/sqrt(2*pi), y);
%sum(w) %Must be close to 1.0
%sum(x.*w) == mean(y)
%sum(x.^(N-2).*w) == mean(y.^(N-2))
%pp = griddedInterpolant(sort(fixedNodes), f(sort(fixedNodes)),'linear');
%sum(pp(x).*w) == mean(pp(x))


function [x,w,nodes] = cappedaccuracy_quad_v3(degree,f,Y)
D = degree(1); N = degree(2);
Kmax = length(Y); %Sample size
nodes = Y(1:D); %Fixed sample set
L = Y(D+1:end); Lsize = length(L);
%Constructing the piecewise linear interpolant
pp = griddedInterpolant(sort(nodes), f(sort(nodes)),'linear');
%-----------Constructing the starting QR-----------
    %Initialize quad rule
    xstart = L(1:N);
    wstart = ones(N,1)/(N);
    %Implicit quad rule
    for ii=N:(Lsize-1)
       %Node addition
       xstart = [xstart;L(ii+1)];
       wstart = [((ii)/(ii+1))*wstart;1./(ii+1)];
       %Update weights
       %V = fliplr(vander(x))';
       V = general_vandermonde(xstart, @(x,k) x.^(k-1), 1:N);
       V(end,:) = pp(xstart);
       nullVec = null(V(1:end,:));
       c = nullVec(:,1);

       % Apply Caratheodory's and determine the node that will be removed
        [alpha, k] = min(wstart(c>0)./c(c>0));
        %[alpha, k] = max(w(c<0)./c(c<0));
        id = find(c > 0);%find(c < 0)
        k = id(k);

        %Node removal
        wstart = wstart-alpha*c;
        xstart(k, :) = []; wstart(k) = [];
    end
%--------Node Embedding-----------------
    iter = Lsize;
    
    %Initialize QR
    x = xstart; w = wstart;
    while iter <= (Kmax)-1
        x(end+1) = nodes(iter+1-Lsize);
        w = [(iter/(iter+1))*w;1/(iter+1)];
        %Construct Vandermonde matrix
        V = general_vandermonde(x, @(x,k) x.^(k-1), 1:N); 
        V(end,:) = pp(x);
        %Applying Caratheodory
        nullVec = null(V);
        c = nullVec(:,1);
%         if size(nullVec,2) > 1
%             size(nullVec,2)
%         end
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




end
                                                                                                                                                                            
