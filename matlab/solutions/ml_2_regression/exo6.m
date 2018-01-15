clf; hold on;
bar(1:p, abs(wSparse), 'b');
bar(1:p, -abs(wRidge), 'r');
legend('Sparse', 'Ridge');
axis tight; box on;
SetAR(1/2);
