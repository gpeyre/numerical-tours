f1 = fw;
clf;
for j=Jmin:Jmax
    a = f1(1:2^j);
    d = f1(2^j+1:2^(j+1));
    a = cconvol(upsampling(a,1),reverse(h),1);
    d = cconvol(upsampling(d,1),reverse(g),1);
    f1(1:2^(j+1)) = a + d;
    j1 = Jmax-j;
    if j1<4
        subplot(4,1,j1+1);
        plot(1:2^j1:N,f1(1:2^(j+1)),'.-'); axis('tight');
        title(['Partial reconstruction, j=' num2str(j)]);
    end
end
