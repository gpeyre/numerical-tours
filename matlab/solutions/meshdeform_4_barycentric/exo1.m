C = zeros(n,n,k);
for i=1:k
    vi = V(:,i);
    U = repmat(vi,[1 n^2])-W;
    nb = normalize( U );
    % length
    d = sqrt( sum(U.^2) );
    s = 1;
    for j=mod([i-2,i],k)+1
        vj = V(:,j);
        na = normalize( repmat(vj,[1 n^2])-W );
        % sign
        si = s*sign(crossp(na,nb));
        % angle
        dp = dotp(na,nb);
        theta = si .* acos(clamp(dp,-1,1));
        % add tangent of half angle
        C(:,:,i) = C(:,:,i) + reshape( tan(theta/2) ./ d, [n n]);
        s = -s;
    end
end
