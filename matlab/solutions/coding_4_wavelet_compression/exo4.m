h = h+1e-10; h = h/sum(h);
E = -sum( h.*log2(h) );
h0 = h0+1e-10; h0 = h0/sum(h0);
E0 = -sum( h0.*log2(h0) );
disp(['Pixels entropy:  ' num2str(E0,2)]);
disp(['Wavelet entropy: ' num2str(E,2)]);
