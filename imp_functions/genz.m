%Function for the six types of multi-dimensional Genz functions with their exact integral over the range [-1,1]^d



function [fn,exact] = genz(x,a,u,type)
d = size(x,2);
npts = size(x,1);
switch type
    case 1
        %Oscillatory
        fn = cos(2*pi*u(1) + x*a);
        exact = (2^d/prod(a))*cos(0.5*sum(a) + 2*pi*u(1))*prod(sin(0.5*a));
    case 2
        %product peak
        unew = repmat(u',[npts,1]);
        anew = repmat(a',[npts,1]);
        fn = prod(1./((anew.^(-2) + (x - unew).^2)),2);
        exact = prod(a.*(atan(a.*u) + atan(a - a.*u)));
    case 3
        %Corner peak
        fn = (1 + x*a).^(-d-1);
        combos = [0,1];
        for ii=1:(d-1)
            combos = combvec([0,1],combos);
        end
        combos = combos';
        exact = prod(1./(1:d))*(1./prod(a))*sum(((-1).^(sum(combos,2))).*(1 + combos*a).^(-1));
    case 4
        %Gaussian
        unew = repmat(u',[npts,1]);
        fn = exp(-((x - unew).^2)*(a).^2);
        exact = ((sqrt(pi)/2)^d)*prod(erf(a.*(1-u)) + erf(a.*u))/prod(a);
    case 5
        %C0 fn
        ploc = 0.5;
        fn = exp(-(abs(x - ploc)*a));
        exact = (1/prod(a))*prod((1 - exp(-a*ploc)) + (1 - exp(-a*(1-ploc))));
    case 6
        %Discontinuous function
        unew = repmat(u',[npts,1]);
        fn = (prod(x <= unew,2)>0).*exp(x*a);
        exact = prod(1./a)*prod(exp(a.*u) - 1);
        
end
end
