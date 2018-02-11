vmlist = [2 5];
j = 3;
k = 0;
clf;
for i=1:length(vmlist)
    vm = vmlist(i);
    options.h = compute_wavelet_filter('Daubechies',vm*2);
    for s=1:3
        k = k+1;
        fW = zeros(n,n);
        if s==1
            fW(1 + 2^j, 1) = 1;
            str = 'Hor';
        elseif s==2
            fW(1,1 + 2^j) = 1;
            str = 'Vert';
        elseif s==3
            fW(1 + 2^j,1 + 2^j) = 1;
            str = 'Diag';
        end
        f1 = perform_wavortho_transf(fW,1,-1,options);
        [tmp,u] = max(f1(:)); [u,v] = ind2sub([n n], u);
        f1 = circshift(f1,[n/2-u n/2-v]);
        subplot(2,3,k);
        surf(f1);
        shading('interp'); view(3);
        axis('tight'); axis('off');
        lighting('phong');
        camlight;
        colormap(jet(256));
        title( strcat([num2str(vm) ' VM, ' str]) );
    end
end
%EXO