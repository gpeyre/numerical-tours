% bounding boxes
B = max(max(abs(z(:,1:2))));
q = 200;
r = linspace(-B,B,q);
[V,U] = meshgrid(r,r);
z1 = [U(:),V(:)];
% test for different R
Rlist = [1 5 10 40];
clf;
for i=1:length(Rlist)
    R=Rlist(i);
    %
    D = distmat(z(:,1:2),z1);
    [Ds,I] = sort(D,1);
    Y = y(I);
    %
    if R==1
        C = Y(1,:);
    else
        h = hist(Y(1:R,:), 1:k);
        [~,C] = max(h);
    end
    C = reshape(C, [q q]);
    % display
    subplot(2,2,i);
    hold on;
    imagesc(r,r,C');
    for i=1:k
        I = find(y==i);
        plot(z(I,1), z(I,2), '.', 'Color', col{i}, 'MarkerSize', ms);
    end
    axis tight; axis equal; axis off;
    title(['R=' num2str(R)]);
    SetAR(1);
end
colormap jet(256);
