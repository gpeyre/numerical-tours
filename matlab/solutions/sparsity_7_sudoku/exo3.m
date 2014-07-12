Acol = [];
for i=1:n
    for k=1:n
        Z = zeros(n,n,n);
        Z(:,i,k) = 1;
        Acol(end+1,:) = Z(:)';
    end
end
clf;
imageplot(Acol);
