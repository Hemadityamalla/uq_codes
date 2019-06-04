%Function to 'marginalize' a given multidimensional function.
function fmarg = marginalize(fevals, evalpts, doi)

    [C,ia,ic] = unique(evalpts(:,doi),'stable');
    fsum = accumarray(ic,fevals);
    fmarg = fsum/histc(ic,1);
end