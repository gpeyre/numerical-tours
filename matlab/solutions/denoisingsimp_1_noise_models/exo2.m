ntrial = 5000;
nlist = 2.^(5:13);
gauss_max = []; gauss_dev = [];
for i=1:length(nlist)
    w = randn(nlist(i), ntrial);
    m = compute_max( abs(w), 1 );
    gauss_max(i) = mean(m);
    gauss_dev(i) = std(m);
end
gauss_max_th = sqrt(2*log(nlist));
clf;
subplot(2,1,1);
hold on;
plot( gauss_max_th, gauss_max, 'r.-' );
plot( gauss_max_th, gauss_max_th, 'b:');
hold off;
axis('tight');
set_label('sqrt(2*log(n))', 'empirical max');
title('Maximum of n dimensional Gaussian noise');
subplot(2,1,2);
plot( gauss_max_th, gauss_dev, '.-' );
axis('tight');
set_label('sqrt(2*log(n))', 'Deviation from max');
