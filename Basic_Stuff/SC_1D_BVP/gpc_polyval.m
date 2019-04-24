function fvals = gpc_polyval(coeffs, eval_pts)
    fvals = 0.0;
    for ii=1:length(coeffs)
       fvals = fvals + hermite(eval_pts, ii-1)*coeffs(ii,1); 
    end
    
end