err = [];
for i=1:size(fList,3)
    % 
    fW = perform_wavortho_transf(fList(:,:,i),Jmin,+1, options);
    e = compute_l2err(fW,1);
    err = [err e];
end
sel = 10:round(n*n/10);
Err = err(sel,:);
Err = Err ./ repmat(Err(1,:), [size(Err,1), 1]);
clf;
h = plot(log10(sel), log10(Err));
if using_matlab()  
    set(h, 'LineWidth', lw);
end
title('log_{10}( \epsilon^2[M]/||f||^2 )');
legend('regular', 'phantom', 'lena', 'mandrill');
axis([log10(sel(1)) log10(sel(length(sel))) -7 0]);
if using_matlab()  
    box('on');
end
