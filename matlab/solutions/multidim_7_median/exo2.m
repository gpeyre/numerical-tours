med = mean(X, 2);
niter = 8; energy = [];
for i=1:niter
    % comute the distance from med to the points
    dist = sqrt( sum( (X-repmat(med,[1 w*w])).^2 ) );
    % compute the weight, take care of not dividing by 0
    weight = 1./max( dist, 1e-10 ); weight = weight/sum(weight);
    % compute the weighted mean
    med = sum( repmat(weight,[3 1]).*X, 2 );
    energy(end+1) = sum( dist );
end
