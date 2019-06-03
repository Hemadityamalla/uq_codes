function [X,W] = multiD_adaptive(f, N, D, S, Y)
    Kmax = length(Y);
    k_level=1;
    nodes = Y(1:D); %Fixed sample set
    L = setdiff(Y,nodes,'stable'); Lsize = length(L);
    %Constructing the piecewise linear interpolant
    pp = griddedInterpolant(sort(nodes), f(sort(nodes)),'linear');


end