vmlist = 2:5;
j = 3;
clf;
for i=1:length(vmlist)
    vm = vmlist(i);
    options.h = compute_wavelet_filter('Daubechies',vm*2);
    fw = zeros(N,1);
    fw(1 + 2^j) = 1;
    f = perform_wavortho_transf(fw,1,-1,options);
    f = circshift(f,N/2);
    subplot(length(vmlist),1,i);
    hh = plot(f); axis('tight');
    if using_matlab()
        set_linewidth(hh, 2);
    end
    title([num2str(vm) ' VM']);
end
%EXO