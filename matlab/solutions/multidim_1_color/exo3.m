pca1rgb = @(g,V,m)applymat(g,V') + repmat(m, [n n 1]);
f1 = pca1rgb(g1,V,m);
% display
clf;
imageplot(f, 'f', 1,2,1);
imageplot(clamp(f1), 'f_1', 1,2,2);
