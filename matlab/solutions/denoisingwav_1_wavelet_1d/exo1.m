Tlist = linspace(.5,6,80)'*sigma;
snr0 = arrayfun( @(T)snr(f0, Theta0W(f,T)), Tlist );
snr1 = arrayfun( @(T)snr(f0, Theta1W(f,T)), Tlist );
clf; plot(Tlist/sigma, [snr0 snr1], 'LineWidth', 2);
axis tight; set_label('T/\sigma', 'SNR');
legend('Hard', 'Soft');
% retrieve best
[~,i] = max(snr0); fBest0 = Theta0W(f,Tlist(i));
[~,i] = max(snr1); fBest1 = Theta1W(f,Tlist(i));
