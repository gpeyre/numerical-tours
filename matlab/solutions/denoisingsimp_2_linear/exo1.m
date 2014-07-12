mu_list = linspace(.5,6,6);
clf;
for i=1:length(mu_list)
    mu = mu_list(i);
    subplot(3,2,i);
    plot( denoise(y,mu) );
    title(['\mu=' num2str(mu)]);
    axis([1 N -.05 1.05]);
end
