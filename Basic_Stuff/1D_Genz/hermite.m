function herm = hermite(pts, degree)
%Generate the probabilist hermite polynomials at the given 'pt' points and
%of a given 'degree'. Here the three term recursion is used

%Make sure that pts is a column vector!!!!!!!

if nargin == 1
    degree = 1; %Gives a default first degree hermite polynomial
end



if degree == -1
    herm = zeros(length(pts),1);
elseif degree == 0
    herm = ones(length(pts),1);
else
    herm = pts.*hermite(pts,degree-1) - (degree-1)*hermite(pts,degree-2);
end


end