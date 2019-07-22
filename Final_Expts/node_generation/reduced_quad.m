function [xred,wred] = reduced_quad(x,w)
N = length(w);

%Constructing the vandermonde matrix excluding the last row
V = vander(x);
V = V';

%Constructing the null vector
 nullVec = null(V(1:end-1,:));
 c = nullVec(:,1);

% Apply Caratheodory's and determine the node that will be removed
[alpha, k] = min(w(c>0)./c(c>0));
id = find(c > 0);
k = id(k);
%Node removal
wred = w-alpha*c;
xred = x;
xred(k, :) = []; wred(k) = [];

end