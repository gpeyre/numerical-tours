Jmin = 0;
fw = f;
clf;
for j=J:-1:Jmin
    fw(1:2^(j+1),1:2^(j+1)) = haar(fw(1:2^(j+1),1:2^(j+1)));
    %
    j1 = J-j;
    if j1<4
        A = fw(1:2^(j+1),1:2^(j+1));
        imageplot(A(1:2^j,2^j+1:2^(j+1)), ['Horizontal, j=' num2str(j)], 3,4, j1 + 1);
        imageplot(A(2^j+1:2^(j+1),1:2^j), ['Vertical, j=' num2str(j)], 3,4, j1 + 5);
        imageplot(A(2^j+1:2^(j+1),2^j+1:2^(j+1)), ['Diagonal, j=' num2str(j)], 3,4, j1 + 9);
    end
end
