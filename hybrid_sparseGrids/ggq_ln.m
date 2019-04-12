function [x,w] = ggq_ln(n)
%This function can be vastly sped up. You need to load and generate rules
%everytime.....but this seems better than having multiple data files.
quad = csvread('ggq_ln_40.csv'); N = 40;
f = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2));
x = quad(:,1);w = quad(:,2);
for ii=N:-1:(n+1)
    [x,w] = reduced_quad(x,w,f); 
end

end