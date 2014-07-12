% compute the threshold T
f = load_signal(name, N);
m_list = [25 50 100]; % number of kept coefficients
options.h = compute_wavelet_filter('Daubechies',4);
clf;
for i=1:length(m_list)
    m = m_list(i);
    fw = perform_wavortho_transf(f,1,+1,options);
    % select threshold
    v = sort(abs(fw(:)));
    if v(1)<v(N)
        v = reverse(v);
    end
    T = v(m);
    fwT = fw .* (abs(fw)>=T);
    % inverse
    f1 = perform_wavortho_transf(fwT,1,-1,options);
    % display
    subplot(length(m_list),1,i);
    hh = plot(f1); axis('tight');
    if using_matlab()
        set_linewidth(hh, 2);
    end
    title(['M=' num2str(m) ', SNR=' num2str(snr(f,f1),3) 'dB']);
end
