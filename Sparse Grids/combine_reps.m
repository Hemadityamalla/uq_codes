function [xnew,wnew] = combine_reps(x,w)
 %Combining repeated entries
 snodes = x';
 sweights = w';
 %1) Homogenizing the zero entries
 snodes(abs(snodes) < 1e-16) = 0;
 %2) Obtaining index list of non-unique nodes
 [C,ia,ic] = unique(snodes,'rows','stable');
 %3) Obtaining non-unique indices(NOT MEMORY EFFICIENT)
%  temp = repmat(ic,[1,length(ic)]);
%  temp2 = tril(temp==temp',-1);
%  [~,~,non_uq_idx] = find(temp(temp2));
%   non_uq_idx = unique(non_uq_idx,'stable');
  %------MEMORY EFFICIENT(at least better than the above implementation)
  non_uq_idx = find(accumarray(ic,1)-1);
 %4) Deleting the non-unique nodes and clubbing the weights
 rep_cache = [];
 for ii=1:length(non_uq_idx)
    reps = find(ic == non_uq_idx(ii));
    rep_cache = [rep_cache;reps(2:end)];
    sweights(reps(1)) = sweights(reps(1)) + sum(sweights(reps(2:end)));
 end
 sweights(rep_cache) = [];
 snodes(rep_cache,:) = [];
 %5)New nodal set and weights
 xnew = snodes';
 wnew = sweights';

end