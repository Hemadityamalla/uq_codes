function [x,w] = smolyakSparseGrid(d,k,growth,quadrature)



N = d+k;%Note: N = d + k, where k is the order of exactness

%Obtaining the combinations of multi-indices
temp = monomialDegrees(d,N);
    %Filtering out the combinations with zero degrees
 temp = temp(find(prod(temp,2)),:);
 
 %Selecting the co1mbinations satisfying N-d+1 <= sum(temp) <= N in the W&W
 %formula
 alpha = temp(sum(temp,2) >= N-d+1 & sum(temp,2)<= N,:);
 
 %For each combination in alpha
 snodes = [];
 sweights = [];
 for aidx=1:size(alpha,1)
     %Below initialization of structs is inefficient as the size of the
     %elements keeps changing
     quadRule(size(alpha,2)).nodes = zeros(1,1);
     quadRule(size(alpha,2)).weights = zeros(1,1);
     %Generating quadrature rules for each combination(how can this be vectorized?)
     [xi,w] = quadgen(growth(alpha(aidx,1)),quadrature);
     quadRule(1).nodes = xi';
     quadRule(1).weights = w';
     %The below initialization is inefficient as the variables are
     %converted to an array whose size keeps increasing.
     tempNodes = quadRule(1).nodes;
     tempWeights = quadRule(1).weights;
     for ii=2:size(alpha,2)
         [xi,w] = quadgen(growth(alpha(aidx,ii)),quadrature);
         quadRule(ii).nodes = xi';
         quadRule(ii).weights = w';
         tempNodes = combvec(tempNodes, quadRule(ii).nodes);
         tempWeights = combvec(tempWeights, quadRule(ii).weights);
     end
     %[tempNodes, tempWeights] = combine_reps(tempNodes, tempWeights);
     %Computing the collective weight for each point (product of the 1D weights and the coeffecient of W&W)
     w  = (-1)^(N - norm(alpha(aidx,:),1)) *nchoosek(d-1,N-norm(alpha(aidx,:),1))*prod(tempWeights,1);
     %Appending the nodes and the weights
     snodes = [snodes,tempNodes];
     sweights = [sweights,w];
 end
%Merging the duplicates
[x, w] = combine_reps(snodes, sweights);

end
