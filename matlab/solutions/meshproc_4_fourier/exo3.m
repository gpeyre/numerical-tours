% nonlinear
v = sort(abs(pvertex(:)).^2, 'ascend');
errN = reverse( cumsum(v) ); 
% linear
v = reverse( abs(pvertex(:)).^2 );
errL = reverse( cumsum(v) ); 
% normalize
errL = errL/errL(1);
errN = errN/errN(1);
% display
close; clf;
plot( log10(1:3*n), log10([errL errN]) ); axis('tight');
axis([0 log10(2*n) -5 0]);
legend('Linear', 'Non-linear');
