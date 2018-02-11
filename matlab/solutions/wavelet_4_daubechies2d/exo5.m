vmlist = 1:3;
clf;
for i=1:length(vmlist)
    vm = vmlist(i);
    options.h = compute_wavelet_filter('Daubechies',vm*2);
    fW = perform_wavortho_transf(f,1,+1,options);
    % select threshold
    v = sort(abs(fW(:)));
    if v(1)<v(n)
        v = reverse(v);
    end
    T = v(m);
    fWT = fW .* (abs(fW)>=T);
    % inverse
    f1 = perform_wavortho_transf(fWT,1,-1,options);
    % display
    imageplot(clamp(f1), strcat([num2str(vm) ' VM, SNR=' num2str(snr(f,f1),3) 'dB']), 1,3,i);
end
