% compute the threshold T
m = 100; % number of kept coefficients
vmlist = 1:3;
clf;
for i=1:length(vmlist)
    vm = vmlist(i);
    options.h = compute_wavelet_filter('Daubechies',vm*2);
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
    subplot(length(vmlist),1,i);
    hh = plot(f1); axis('tight');
    if using_matlab()
        set_linewidth(hh, 2);
    end
    title([num2str(vm) ' VM, SNR=' num2str(snr(f,f1),3) 'dB']);
end
