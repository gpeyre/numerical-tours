mu_list = linspace(.5,6,6);
clf;
for i=1:length(mu_list)
    mu = mu_list(i);
    subplot(2,3,i);
    imageplot( denoise(y,mu) );
    title(['\mu=' num2str(mu)]);
end
