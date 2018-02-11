J = log2(N)-1;
selj = ( J-2:J )-3;
pos = [0 .5];
f = [];
clf;
k = 0;
for j=selj
    k = k+1;
    for q=1:length(pos)
        fw = zeros(N,1);
        p = 1 + (1+pos(q))*2^j;
        fw(p) = 1;
        f(:,q) = perform_haar_transf(fw,1,-1);
        f(:,q) = circshift(f(:,q),N/4);
    end
    f(1:N/2-1,2) = nan(); f(N/2+1:N,1) = nan();
    subplot(3,1,k);
    hh = plot(f); axis('tight');
    axis([1 N min(f(:))*1.05 max(f(:))*1.05]);
    if using_matlab()
        set_linewidth(hh,2);
    end
end
