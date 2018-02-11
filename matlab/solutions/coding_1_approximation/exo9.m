f1 = fL;
for i=1:n/w
    for j=1:n/w
        seli = (i-1)*w+1:i*w;
        selj = (j-1)*w+1:j*w;
        f1(seli,selj) = idct2( f1(seli,selj) );
    end
end
disp(strcat((['Error |f-f1|/|f| = ' num2str(norm(f(:)-f1(:))/norm(f(:)))])));
