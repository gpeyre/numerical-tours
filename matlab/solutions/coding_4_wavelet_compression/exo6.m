Etot = compute_entropy([MWH MWV MWD]);
Ehor = compute_entropy(MWH);
Ever = compute_entropy(MWV);
Edia = compute_entropy(MWD);
disp([ 'Entropy, all:  ' num2str(Etot,3) ]);
disp([ 'Entropy, hor:  ' num2str(Ehor,3) ]);
disp([ 'Entropy, vert: ' num2str(Ever,3) ]);
disp([ 'Entropy, diag: ' num2str(Edia,3) ]);
