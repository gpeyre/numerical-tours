% bounding boxes
B = max(max(abs(Z(:,1:2))));
q = 200;
r = linspace(-B,B,q);
[V,U] = meshgrid(r,r);
z1 = [U(:),V(:)];
% test for different R
Rlist = [1 5 10 40];
clf;
for ir=1:length(Rlist)
    R=Rlist(ir);
    %
    D = distmat(Z(:,1:2),z1);
    [Ds,I] = sort(D,1);
    ys = y(I);
    %
    if R==1
        C = ys(1,:);
    else
        h = hist(ys(1:R,:), 1:k);
        [~,C] = max(h);
    end
    C = reshape(C, [q q]);
    % maps class to color
    Cr = zeros(q,q,3);
    for i=1:k
        for a=1:3
            Cr(:,:,a) = Cr(:,:,a) + (C==i)*col(a,i);
        end
    end
    % display
    subplot(2,2,ir);
    hold on;
    imagesc(r,r,permute(Cr,[2 1 3]));
    for i=1:k
        I = find(y==i);
        plot(Z(I,1), Z(I,2), 'o', 'MarkerFaceColor', col(:,i)*.9, 'MarkerSize', 5, 'MarkerEdgeColor', col(:,i)*.5);
    %    plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', ms);
    end
    axis tight; axis equal; axis off;
    title(['R=' num2str(R)]);
    SetAR(1);
end
colormap jet(256);
