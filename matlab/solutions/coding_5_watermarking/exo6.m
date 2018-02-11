pfa = 10.^( linspace(-10,-2) );
ptp = [];
for i=1:length(pfa)
    T = sqrt(2) * sigma0 * erfinv(1-2*pfa(i));
    ptp(i) = sum(c>T)/length(c);
end
clf;
h = plot(log10(pfa), ptp);
set(h, 'LineWidth', 2);
axis('tight');
xlabel('log_{10}(p_{FA})');
ylabel('p_{TP}');
