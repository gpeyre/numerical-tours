k = 1;      % half window
w = 2*k+1;  % window size
Mmed = M;
clf;
for i=1:8
    Mmed = perform_median_filtering(Mmed,k);
    if mod(i,2)==1 & (i<=7)
        imageplot(clamp(Mmed), strcat(['Iter 5' num2str(i)]), 2,2, (i-1)/2+1 );
    end
    err(i) = snr(M0,Mmed);
end
