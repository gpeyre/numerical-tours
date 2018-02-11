jlist = J-(1:4);
fw = perform_haar_transf(f,1,+1);
clf;
for i=1:length(jlist)
    j = jlist(i);
    fw1 = zeros(n); fw1(1:2^j,1:2^j) = fw(1:2^j,1:2^j);  
    f1 = perform_haar_transf(fw1,1,-1);
    % display
    subplot(2,2,i);
    imageplot(f1);
    title( strcat(['j=' num2str(j) ', SNR=' num2str(snr(f,f1),3) 'dB']) );
end
