jlist = J-(1:3);
fw = perform_haar_transf(f,1,+1);
clf;
for i=1:length(jlist)
    j = jlist(i);
    fw1 = fw; fw1(2^j+1:end) = 0;    
    f1 = perform_haar_transf(fw1,1,-1);
    % display
    subplot(length(jlist),1,i);
    hh = plot(f1); axis('tight');
    if using_matlab()
        set_linewidth(hh,2);
    end
    title( strcat(['j=' num2str(j) ', SNR=' num2str(snr(f,f1),3) 'dB']) );
end
