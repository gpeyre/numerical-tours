Aenc = [];
for i=1:n
    for j=1:n
        Z = zeros(n,n,n);
        Z(i,j,:) = 1;
        Aenc(end+1,:) = Z(:)';
    end
end
clf;
imageplot(Aenc);
