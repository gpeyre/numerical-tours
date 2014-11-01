Jmax = log2(n)-1;
Jmin = 0;
fW = f;
clf;
for j=Jmax:-1:Jmin
    A = fW(1:2^(j+1),1:2^(j+1));
    for d=1:2
        Coarse = subsampling(cconvol(A,h,d),d);
        Detail = subsampling(cconvol(A,g,d),d);
        A = cat3(d, Coarse, Detail );
    end
    fW(1:2^(j+1),1:2^(j+1)) = A;
    j1 = Jmax-j;
    if j1<4
        imageplot(A(1:2^j,2^j+1:2^(j+1)), strcat(['Horizontal, j=' num2str(j)]), 3,4, j1 + 1);
        imageplot(A(2^j+1:2^(j+1),1:2^j), strcat(['Vertical, j=' num2str(j)]), 3,4, j1 + 5);
        imageplot(A(2^j+1:2^(j+1),2^j+1:2^(j+1)), strcat(['Diagonal, j=' num2str(j)]), 3,4, j1 + 9);
    end
end
