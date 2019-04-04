function fn = haar(x,j,k)

PhiW  = @(x)  (1*(x <= 0.5));
wavFam = @(x,j,k) 2^(j/2)*PhiW(x*2^j - k);

fn = wavFam(x,j,k);
end