nj = size(face{1},2);
f1 = fw;
for j=1:J-1
    fj = f1(1:4*nj);
    fj = reshape(fj, [nj 4]);
    fj = fj*U';
    f1(1:4*nj) = fj(:);
    nj = nj*4;
end
e = norm(f-f1)/norm(f);
disp(['Bijectivity deviation (should be 0): ' num2str(e)]);
