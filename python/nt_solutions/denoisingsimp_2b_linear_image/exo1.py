mu_list = linspace(.5,6,6);
clf;
for i in arange(0,len(mu_list)):
    mu = mu_list[i];
    subplot(2,3,i+1);
    imageplot( denoise(y,mu) );
    title('$\mu=' + str(mu) + '$');
