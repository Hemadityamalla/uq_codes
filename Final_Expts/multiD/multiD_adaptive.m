function [X,W] = multiD_adaptive(f, N, D, S, Y, dim)
    Kmax = length(Y);
    nodes = Y(1:D); %Fixed sample set
    L = setdiff(Y,nodes,'stable'); Lsize = length(L);
    %Constructing the tensor product of the fixed sample set
    nodes_tensor = setprod(nodes, dim);
    f_val = f(nodes_tensor); %Here we replace by PDE evaluations
    %Constructing a QR in each dimension
    for i_dim = 1:dim
        [x,w] = fdaqr([N,D],marginalize(f_val,nodes_tensor, i_dim), L, nodes);
        QR{i_dim}.nodes = x'; QR{i_dim}.weights = w';
    end
    if dim == 1
        X = x';W = w';
    else
        X = QR{1}.nodes; W = QR{1}.weights;
        for i_dim=2:dim
            X = combvec(X, QR{i_dim}.nodes); W = combvec(W, QR{i_dim}.weights);
        end
    end
    X = X'; W = prod(W',2);
end