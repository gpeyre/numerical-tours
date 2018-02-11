Ainp = [];
for i=1:length(I)
    Z = zeros(n,n,n);
    Z(I(i), J(i), v(i)) = 1; % double( ==k );
    Ainp(end+1,:) = Z(:)';
end
clf;
imageplot(Ainp);
