dList = [.4 .6 .8 1.2];
clf;
for i=1:length(dList)
    delta = dList(i)/fc;
    x0 = [.5-delta .5 .5+delta]';
    %
    w = ones(N,1);
    Gamma = [];
    for k=0:d
        Gamma = [Gamma, diag(w) * Fourier(fc,x0)];
        % derivate the filter
        w = w .* 2i*pi .* (-fc:fc)';
    end
    %
    pV = pinv(Gamma') * [sign(a0); zeros(n,1)];
    etaV = PhiS(fc, u, pV);
    %
    subplot(2,2,i);
    hold on;
    stem(x0, sign(a0), 'k.--', 'MarkerSize', ms, 'LineWidth', lw);
    plot([0 1],  [1 1], 'k--', 'LineWidth', lw);
    plot([0 1], -[1 1], 'k--', 'LineWidth', lw);
    plot(u, etaV, 'b', 'LineWidth', lw);
    axis([0 1 -1.4 1.4]);
    set(gca, 'XTick', [], 'YTick', [-1 1]);
    box on;
    title(['\delta=' num2str(dList(i)) '/f_c']);
end
