function cheb = chebyshev(pts, degree)
%Generate Chebyshev polynomials of first type at the given 'pts' points and
%of a given 'degree'. Here the three term recursion is used

%Make sure that pts is a column vector!!!!!!!

if nargin == 1
    degree = 1; %Gives a default first degree Chebyshev polynomial
end



if degree == 0
    cheb = ones(length(pts),1);
elseif degree == 1
    cheb = pts;
else
    cheb = 2*pts.*chebyshev(pts,degree-1) - chebyshev(pts,degree-2);
end


end