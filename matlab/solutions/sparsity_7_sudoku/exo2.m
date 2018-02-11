Arow = [];
for i=1:n
    for k=1:n
        Z = zeros(n,n,n);
        Z(i,:,k) = 1;
        Arow(end+1,:) = Z(:)';
    end
end
clf;
imageplot(Arow);
