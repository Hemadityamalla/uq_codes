function herm = hermite(pts, degree)
%Generate the probabilist hermite polynomials at the given 'pt' points and
%of a given 'degree'. Here the three term recursion is used
%degree is usually a column of size (number of dimensions)
%Make sure that pts is a matrix of size (Number of points*number of dimensions)
[Npts, d] = size(pts);
%if nargin == 1
%    degree = 1; %Gives a default first degree hermite polynomial
%end
if d == 1
    if degree == -1
        herm = zeros(length(pts),1);
    elseif degree == 0
        herm = ones(length(pts),1);
    else
        herm = pts.*hermite(pts,degree-1) - (degree-1)*hermite(pts,degree-2);
    end
    return;
end


%The code below can be vectorized-- future work
polymat = zeros(Npts, d);
for i_pts = 1:Npts
    for i_d = 1:d
       polymat(i_pts,i_d) = hermite(pts(i_pts,i_d),degree(i_d)); 
    end
end

herm = prod(polymat,2);

end