function leg = legendre(pts, degree)
%Generate the probabilist hermite polynomials at the given 'pt' points and
%of a given 'degree'. Here the three term recursion is used

%Make sure that pts is a column vector!!!!!!!

if nargin == 1
    degree = 1; %Gives a default first degree hermite polynomial
end



if degree == 0
    leg = ones(length(pts),1);
elseif degree == 1
    leg = pts;
else
    leg = ((2*degree - 1)/(degree))*pts.*legendre(pts,degree-1) - ((degree-1)/degree)*legendre(pts,degree-2);
end


end