function [xred,wred] = reduced_quad(x,w,f,selection)
if nargin==3
    selection = 'none';
end
N = length(w);

%coeffs for the ggq of J.Ma et.al. paper--in this case, the order is 2N,
%but we take only upto N
coeffs = 1:N;

%Constructing the vandermonde matrix excluding the last row
V = general_vandermonde(x, f,coeffs(1:end-1));

%Constructing the null vector
[~,~,U] = svd(V);
nnv = diff(size(V));%size(V,2) - size(V,1);
nullVec = U(:,end-nnv:end);
c = nullVec(:,end);
% nullVec = null(V);
% c = nullVec(:,1);

% Apply Caratheodory's and determine the node that will be removed
if strcmp(selection,'none')
[alpha, k] = min(w(c>0)./c(c>0));
%[alpha, k] = max(w(c<0)./c(c<0));
id = find(c > 0);%find(c < 0)
k = id(k);
else
    [alpha1, k1] = min(w(c>0)./c(c>0));
    [alpha2, k2] = max(w(c<0)./c(c<0));
    id1 = find(c > 0);
    id2 = find(c < 0);
    k1 = id1(k1); k2 = id2(k2);
end

%Node removal
wred = w-alpha*c;
xred = x;
xred(k, :) = []; wred(k) = [];

end