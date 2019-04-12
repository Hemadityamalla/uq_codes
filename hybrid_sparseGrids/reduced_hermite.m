function [x,w] = reduced_hermite(n)
%This function can be vastly sped up. You need to load and generate rules
%everytime.....but this seems better than having multiple data files.
f = @(x,k) x.^(k-1);
N=15;
[x,w] = gaussQuad(N,'Hermite');
for ii=N:-1:(n+1)
    [x,w] = reduced_quad(x,w,f); 
end

end