fL = zeros(n,n);
for i=1:n/w
    for j=1:n/w
        seli = (i-1)*w+1:i*w;
        selj = (j-1)*w+1:j*w;
        fL(seli,selj) = dct2( f(seli,selj) );
    end
end
