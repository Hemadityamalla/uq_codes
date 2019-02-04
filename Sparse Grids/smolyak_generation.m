clear;clc; close all;

d = 10; k = 3; %Level of the smolyak quadrature
N = d+k;%Note: N = d + k, where k is the order of exactness
polyBasis = 'Legendre';
%Obtaining the combinations of multi-indices
temp = monomialDegrees(d,N);
    %Filtering out the combinations with zero degrees
 temp = temp(find(prod(temp,2)),:);
 
 %Selecting the combinations satisfying N-d+1 <= sum(temp) <= N in the W&W
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
     [xi,w] = gaussQuad(alpha(aidx,1),polyBasis);
     quadRule(1).nodes = xi';
     quadRule(1).weights = w';
     %The below initialization is inefficient as the variables are
     %converted to an array whose size keeps increasing.
     tempNodes = quadRule(1).nodes;
     tempWeights = quadRule(1).weights;
     for ii=2:size(alpha,2)
         [xi,w] = gaussQuad(alpha(aidx,ii),polyBasis);
         quadRule(ii).nodes = xi';
         quadRule(ii).weights = w';
         tempNodes = combvec(tempNodes, quadRule(ii).nodes);
         tempWeights = combvec(tempWeights, quadRule(ii).weights);
     end
     %Computing the collective weight for each point (product of the 1D weights and the coeffecient of W&W)
     w  = (-1)^(N - norm(alpha(aidx,:),1)) *nchoosek(d-1,N-norm(alpha(aidx,:),1))*prod(tempWeights,1);
     %Appending the nodes and the weights
     snodes = [snodes,tempNodes];
     sweights = [sweights,w];
 end

 %Combining repeated entries
 snodes = snodes';
 sweights = sweights';
 %1) Homogenizing the zero entries
 snodes(abs(snodes) < 1e-16) = 0;
 %2) Obtaining index list of non-unique nodes
 [C,ia,ic] = unique(snodes,'rows','stable');
 %3) Obtaining non-unique indices(NOT MEMORY EFFICIENT)
 temp = repmat(ic,[1,length(ic)]);
 temp2 = tril(temp==temp',-1);
 [~,~,non_uq_idx] = find(temp(temp2));
 %4) Deleting the non-unique nodes and clubbing the weights
 non_uq_idx = unique(non_uq_idx,'stable');
 rep_cache = [];
 for ii=1:length(non_uq_idx)
    reps = find(ic == non_uq_idx(ii));
    rep_cache = [rep_cache;reps];
    sweights(reps(1)) = sweights(reps(1)) + sum(sweights(reps(2:end)));
 end
 sweights(rep_cache) = [];
 snodes(rep_cache,:) = [];
 %5)New nodal set and weights
 snodes = snodes';
 sweights = sweights';



%Expt before exams----------------
%  %Taking care of the repeated nodes (0,0,0...0) due to finite precision
%  %arithmetic
%  snodes = snodes'; sweights = sweights';
%  snodes(abs(snodes) < 1e-15) = 0;
%  [finalquad, ic, ia] = unique(snodes,'rows','stable');
%  ioi = find(ia==find(accumarray(ia,1)>1));
%  sweights(ioi(1)) = sweights(ioi(1)) + sum(sweights(ioi(2:end)));%Adding the weights of the similar entries
%  sweights(ioi(2:end)) = []; %Deleting the individual weights
%  
%  %Now the final quadrature rule is (finalquad, sweights)
%  %snodes = finalquad;
 


 % This part is just for visualization, works only for d=2
if d == 2
    scatter(snodes(1,:),snodes(2,:),'.');
    lim = 1;
    xlim([-lim,lim]);ylim([-lim,lim]);zlim([-lim,lim])
end