Tlist = linspace(.5,2,20)*sigma;
snr_stein = arrayfun( @(T)snr(f0, ThreshWav(f,T)), Tlist );
clf; plot(Tlist/sigma, snr_stein(:), 'LineWidth', 2);
axis tight; set_label('T/\sigma', 'SNR');
% retrieve best
[~,i] = max(snr_stein); T = Tlist(i);
fTI = ThreshWav(f,T);
