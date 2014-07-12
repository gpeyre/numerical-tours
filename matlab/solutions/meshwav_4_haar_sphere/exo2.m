fw = f;
nj = length(f);
for j=1:J-1
    fj = fw(1:nj);
    fj = reshape(fj, [nj/4 4]);
    fj = fj*U;
    fw(1:nj) = fj(:);
    nj = nj/4;
end
% Check for orthogonality.
e = abs(norm(f)-norm(fw))/norm(f);
disp(['Orthogonality deviation (should be 0): ' num2str(e)]);
