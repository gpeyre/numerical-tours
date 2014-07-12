Q = 100;
[h0,t] = hist(f(:), Q);
[h1,t] = hist(g(:), Q);
clf;
subplot(2,1,1);
bar(t,h0*Q/n^2); axis([0 1 0 6]);
subplot(2,1,2);
bar(t,h1*Q/n^2); axis([0 1 0 6]);
