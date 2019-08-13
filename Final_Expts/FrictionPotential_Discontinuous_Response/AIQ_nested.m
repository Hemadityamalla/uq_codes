
function [QR] = AIQ_nested(degree,level,Y)
    D = degree(1); N = degree(2);
    Kmax = length(Y); %Sample size
                                                        %---------First level--------------
    k_level=1;
    nodes = Y(1:D); %Fixed sample set
    L = setdiff(Y,nodes,'stable');
    [x,w] = fdaqr([N,D], Potential_friction_motion(nodes), L, nodes);
    
    %Storing the generated QR
    QR{k_level}.nodes = x; QR{k_level}.weights = w;
                                        %---------------Subsequent levels-------------
    k_level = 2;
    while k_level <= level
        degreeInt = N + (k_level-1);
        nodes = x; %Fixed sample set
        L = setdiff(Y,nodes,'stable');
        [x,w] = fdaqr([degreeInt, length(x)], Potential_friction_motion(nodes), L, nodes);
        QR{k_level}.nodes = x; QR{k_level}.weights = w;
        k_level = k_level+1;
    end



end
                                                                                                                                                                            
