cR = sort(abs(fL(:)));
if cR(n^2)>cR(1)
    cR = reverse(cR); % be sure it is in reverse order
end
err_ldct = max( norm(f(:))^2 - cumsum(cR.^2), 1e-10); 
clf; hold('on');
h = plot(log10(err_fft / norm(f(:))^2), 'r'); 
if using_matlab()  
    set(h, 'LineWidth', lw);
end
h = plot(log10(err_dct / norm(f(:))^2), 'g'); 
if using_matlab()  
    set(h, 'LineWidth', lw);
end
h = plot(log10(err_ldct / norm(f(:))^2), 'b'); 
if using_matlab()  
    set(h, 'LineWidth', lw);
end
h = plot(log10(err_wav / norm(f(:))^2), 'k'); 
if using_matlab()  
    set(h, 'LineWidth', lw);
end
axis(ax);
title('log_{10}( \epsilon^2[M]/||f||^2 )');
legend('Fourier', 'DCT', 'local-DCT', 'Wavelets');
if using_matlab()  
    box('on');
end
