err_fft = max( norm(f(:))^2 - cumsum(cR.^2), 1e-10); 
clf;
h = plot(log10(err_fft / norm(f(:))^2)); 
if using_matlab()  
    set(h, 'LineWidth', lw);
end
ax = [1 n^2/50 -2.35 0];
axis(ax);
title('log_{10}( \epsilon^2[M]/||f||^2 )');
