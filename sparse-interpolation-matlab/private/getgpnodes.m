function x = getgpnodes(allnx)
% Generate Gauss-Patterson nodes for barypdstepgp
%    (internal function)

d = uint8(length(allnx));
x = zeros(sum(allnx),1);

aid = uint32(0);
for k = 1:d
  x(1+aid:aid+allnx(k)) = gpabsc(log2(double(allnx(k))+1)-1);
	aid = aid + allnx(k);
end
