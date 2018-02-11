% init
I = randperm(n);
I = find(y==1);
I = I(1:k);
C = X(I,:);
niter = 16;
it_dist = 0;
clf;
for it=1:niter
    % NN
    D = distmat(X,C);
    [~,yb] = min(D, [], 2);
    % display
    if it<=3 || it==niter
        it_dist = it_dist+1;
        subplot(2,2,it_dist);
        hold on;
        for i=1:k
            I = find(yb==i);
            plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', 25);
        end
        CV = (C-repmat(mean(X,1), [k 1]))*V;
        for i=1:k
            plot(CV(i,1), CV(i,2), 'o', 'MarkerFaceColor', col(:,i), 'MarkerSize', 12, 'MarkerEdgeColor', 'k');
        end
        axis tight; axis equal; axis off; title(['Iter 3' num2str(it)]);
        SetAR(1);
        % update centroids
        for l=1:k
            C(l,:) = mean( X(yb==l,:), 1 );
        end
    end
end
