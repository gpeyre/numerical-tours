Jmax = log2(N)-1; Jmin = 0;
fw = f;
clf;
subplot(4,1,1);
plot(f); axis('tight'); title('Signal');
for j=Jmax:-1:Jmin
    a1 = fw(1:2^(j+1));
    a = subsampling(cconvol(a1,h));
    d = subsampling(cconvol(a1,g));
    fw(1:2^(j+1)) = cat(1, a, d );
    j1 = Jmax-j;
    if j1<3
        subplot(4,1,j1+2);
        plot(1:2^(j1+1):N,d);  axis('tight');
        title(['Details, j=' num2str(j)]);
    end    
end
