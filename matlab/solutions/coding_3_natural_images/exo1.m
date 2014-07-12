[v1,I1] = sort(M1(:));
[v2,I2] = sort(M2(:));
Meq = M1;
Meq(I1) = v2;
Meq = reshape(Meq, [n n]);
clf;
imageplot(M1, 'Original', 1,2,1);
imageplot(Meq, 'Equalized', 1,2,2);
