% compute the threshold T
m_list = round([.05 .1 .2]*N); % number of kept coefficients
fw = perform_haar_transf(f,1,+1);
clf;
for i=1:length(m_list)
    m = m_list(i);
    % select threshold
    v = sort(abs(fw(:)));
    if v(1)<v(N)
        v = reverse(v);
    end
    T = v(m);
    fwT = fw .* (abs(fw)>=T);
    % inverse
    f1 = perform_haar_transf(fwT,1,-1);
    % display
    subplot(length(m_list),1,i);
    hh = plot(f1); axis('tight');
    if using_matlab()
        set_linewidth(hh,2);
    end
    title( strcat(['m=' num2str(m) ', SNR=' num2str(snr(f,f1),3) 'dB']) );
end
