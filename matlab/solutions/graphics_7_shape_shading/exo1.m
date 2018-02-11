v = 2*rand(3,4)-1; v(3,:) = abs(v(3,:));
v = v ./ repmat( sqrt(sum(v.^2)), [3 1] );
clf;
for i=1:4
    d = v(:,i);
    L = sum( N .* repmat(reshape(d,[1 1 3]), [n n 1]),3 );
    subplot(2,2,i);
    imageplot(max(L,vmin));
end
