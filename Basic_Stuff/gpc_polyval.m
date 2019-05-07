function fvals = gpc_polyval(coeffs, eval_pts,polytype)
    fvals = 0.0;
    if strcmp(polytype, 'hermite')
    for ii=1:length(coeffs)
       fvals = fvals + hermite(eval_pts, ii-1)*coeffs(ii,1); 
    end
    elseif strcmp(polytype, 'legendre')
       for ii=1:length(coeffs)
       fvals = fvals + legendre(eval_pts, ii-1)*coeffs(ii,1); 
       end
    else
        error('wrong polytype in polyval');
    end
    
end