function [X,W] = multiD_adaptive(f, N, D, S, Y, dim)
    Kmax = length(Y);
    k_level=1;
    nodes = Y(1:D); %Fixed sample set
    L = setdiff(Y,nodes,'stable'); Lsize = length(L);
    %Constructing the tensor product of the fixed sample set
    nodes_tensor = setprod(nodes, dim);
    f_val = f(nodes_tensor);
    %Constructing a QR in each dimension
    for i_dim = 1:dim
        
    end
    pp = griddedInterpolant(sort(nodes), f(sort(nodes)),'linear');


end