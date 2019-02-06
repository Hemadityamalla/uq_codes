function [x,w] = quadgen(n,q)

if strcmp(q, 'ClenshawCurtis')
    [x,w] = clencurt(n);
else
    [x,w] = gaussQuad(n,q);
end

end