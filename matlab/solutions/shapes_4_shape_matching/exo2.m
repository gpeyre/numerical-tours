C = zeros(nbound);
for i=1:nbound
    for j=1:nbound
        C(i,j) = norm(D{1}(:,i)-D{2}(:,j));
    end
end
clf;
imageplot(C);
