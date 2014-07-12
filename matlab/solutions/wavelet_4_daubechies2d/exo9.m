Jmax = log2(n)-1;
options.h = compute_wavelet_filter('Daubechies',6);
selj = ( Jmax-2:Jmax )-3;
clf;
k = 0;
for j=selj
    for s=1:3
        k = k+1;
        fW = zeros(n,n);
        MWj = zeros(2^j,2^j);
        MWj(1,1) = 1;
        if s==1
            fW(1:2^j,2^j+1:2*2^j) = MWj;
            str = 'Hor';
        elseif s==2
            fW(2^j+1:2*2^j,1:2^j) = MWj;
            str = 'Vert';
        elseif s==3
            fW(2^j+1:2*2^j,2^j+1:2*2^j) = MWj;
            str = 'Diag';
        end
        f1 = perform_wavortho_transf(fW,1,-1,options);
        f1 = circshift(f1,[n n]/2);
        imageplot(f1, strcat(['j=' num2str(j) ' ' str]), 3,3,k);
    end
end
