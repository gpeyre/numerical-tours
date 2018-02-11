k = 100;
t = linspace(0,.5,k);
clf;
subplot(2,1,1);
bar(t, histc(f0(:), t)*k/n^2); axis tight;
title('Input');
subplot(2,1,2);
bar(t, histc(f(:), t)*k/n^2);  axis tight;
title('Synthesized');
