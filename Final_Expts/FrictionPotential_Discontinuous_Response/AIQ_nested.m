
function [QR] = AIQ_nested(degree,level,Y)
    D = degree(1); N = degree(2);
    Kmax = length(Y); %Sample size
                                                      %---------First level--------------
    k_level=1;
    nodes = Y(1:D); %Fixed sample set
    L = setdiff(Y,nodes,'stable');
    QR{k_level}.fnevals = Potential_friction_motion(nodes);
    [x,w] = fdaqr([N,D], QR{k_level}.fnevals, L, nodes);
    %Storing the generated QR 
    QR{k_level}.nodes = x; QR{k_level}.weights = w; 
                                        %---------------Subsequent levels-------------
    k_level = 2;
    while k_level <= level
        degreeInt = N + (k_level-1);
        [newnodes, idx] = setdiff(x, nodes,'stable');
        leftovers = setdiff((1:length(x))',idx,'stable');
        QR{k_level}.fnevals = zeros(length(x),1);
        newfnevals = Potential_friction_motion(newnodes);
        QR{k_level}.fnevals(idx) = newfnevals;
        QR{k_level}.fnevals(leftovers) = QR{k_level-1}.fnevals;
        nodes = x; %Fixed sample set
        L = setdiff(Y,nodes,'stable');
        [x,w] = fdaqr([degreeInt, length(x)], QR{k_level}.fnevals, L, nodes);
        QR{k_level}.nodes = x; QR{k_level}.weights = w;QR{k_level-1}.QoI = QR{k_level}.fnevals;
        k_level = k_level+1;
    end
    %Generating QoI array for the last level
    [newnodes, idx] = setdiff(x, nodes,'stable');
    leftovers = setdiff((1:length(x))',idx,'stable');
    newfnevals = Potential_friction_motion(newnodes);
    QR{k_level-1}.QoI = zeros(length(x),1);
    QR{k_level-1}.QoI(idx) = newfnevals;
    QR{k_level-1}.QoI(leftovers) = QR{k_level-1}.fnevals;
end
                                                                                                                                                                            
