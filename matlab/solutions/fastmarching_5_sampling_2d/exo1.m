% compute quad of values
V = [];
v = Q(1:end-1,1:end-1); V = [V v(:)];
v = Q(2:end,1:end-1); V = [V v(:)];
v = Q(1:end-1,2:end); V = [V v(:)];
v = Q(2:end,2:end); V = [V v(:)];
V = sort(V,2);
V = unique(V, 'rows');
d = (V(:,1)~=V(:,2)) + (V(:,2)~=V(:,3)) + (V(:,3)~=V(:,4));
V = V(d==2,:);
for i=1:size(V,1)
    V(i,1:3) = unique(V(i,:));
end
faces = V(:,1:3)';
