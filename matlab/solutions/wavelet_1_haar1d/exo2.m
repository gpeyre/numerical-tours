f1 = fw;
clf;
for j=Jmin:J
    f1(1:2^(j+1)) = ihaar(f1(1:2^j), f1(2^j+1:2^(j+1)));
    j1 = J-j;
    if j1<4
        subplot(4,1,j1+1);
        plot(1:2^j1:N,f1(1:2^(j+1)),'.-'); axis('tight');
        title( strcat(['Partial reconstruction, j=' num2str(j)]) );
    end
end
