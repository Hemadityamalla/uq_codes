function [x,w] = quadgen(n,q)

if strcmp(q, 'ClenshawCurtis')
    [x,w] = clencurt(n);
elseif strcmp(q,'Legendre')
    [x,w] = gaussQuad(n,q);
elseif strcmp(q,'RHermite')
    [x,w] = reduced_hermite(n);
elseif strcmp(q,'RLegendre')
    [x,w] = reduced_legendre(n);
elseif strcmp(q,'RClenshaw')
    [x,w] = reduced_cc(n);
else
    if n <= 40
        [x,w] = ggq_ln(n);
    else
        error('Cant go beyond 40 nodes....yet.');
    end

end