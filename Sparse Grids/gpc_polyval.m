function fvals = gpc_polyval(coeffs, eval_pts, type, ordering)
    
    fvals = 0.0;
    if nargin == 3 %univariate case
        if strcmp(type,'hermite')
            for ii=1:length(coeffs)
               fvals = fvals + hermite(eval_pts, ii-1)*coeffs(ii,1); 
            end
        elseif strcmp(type,'legendre')
            for ii=1:length(coeffs)
           fvals = fvals + legendre(eval_pts, ii-1)*coeffs(ii,1); 
            end
        else
            error('Polynomial type not given! OR Wrong type given');
        end
    elseif nargin == 4 %Multivariate case
        if strcmp(type,'hermite')
            for ii=1:length(coeffs)
               fvals = fvals + hermite(eval_pts, ordering(ii,:)')*coeffs(ii,1); 
            end
        elseif strcmp(type, 'legendre')
            for ii=1:length(coeffs)
               fvals = fvals + legendre(eval_pts, ordering(ii,:)')*coeffs(ii,1); 
            end
        else
            error('Polynomial type not given! OR Wrong type given');
        end
    
            
    end
    
end