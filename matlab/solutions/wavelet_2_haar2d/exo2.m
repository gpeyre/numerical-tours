f1 = fw;
clf;
for j=Jmin:J
    s = 1:2^j; t = 2^j+1:2^(j+1); u = 1:2^(j+1);
    f1(u,u) = ihaar(f1(s,s),f1(s,t),f1(t,s),f1(t,t)); 
    %
    j1 = J-j;
    if j1>0 & j1<5
        A = f1(1:2^(j+1),1:2^(j+1));
        subplot(2,2,j1);
        imageplot(A, ['Partial reconstruction, j=' num2str(j)]);
    end
end
