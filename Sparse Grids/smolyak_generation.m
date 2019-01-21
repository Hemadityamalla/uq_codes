clear;clc;

d = 2; N = 6;
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
 
 % This part is just for visualization, works only for d=2
if d == 2
    scatter(snodes(1,:),snodes(2,:),'+');
    lim = 1;
    xlim([-lim,lim]);ylim([-lim,lim]);zlim([-lim,lim])
end