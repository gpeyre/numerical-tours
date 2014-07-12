Mmed = zeros(n,n,3);
niter = 15; t = 0;
for x=1:n
    selx = x-k:x+k; selx = mod(selx-1,n)+1;
    for y=1:n
        t = t+1; % progressbar(t,n*n);
        % extract patch
        sely = y-k:y+k; sely = mod(sely-1,n)+1;
        X = reshape( M(selx,sely,:), [w*w 3])';
        % compute median
        med = mean(X, 2);
        for i=1:niter
            dist = sqrt( sum( (X-repmat(med,[1 w*w])).^2 ) );
            weight = 1./max( dist, 1e-10 ); weight = weight/sum(weight);
            med = sum( repmat(weight,[3 1]).*X, 2 );
        end
        % store result
        Mmed(x,y,:) = reshape(med, [1 1 3]);
    end
end
pmed = snr(M0,Mmed);
% display the results
clf;
imageplot(clamp(Mindep), strcat(['1D median, SNR=' num2str(pindep)]), 1,2,1 );
imageplot(clamp(Mmed), strcat(['3D median, SNR=' num2str(pmed)]), 1,2,2 );
