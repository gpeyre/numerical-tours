F = zeros(n,n,2);
for i=1:m
    for j=1:m
        % locations
        x = (i-1)*w+1;
        y = (j-1)*w+1;
        selx = clamp( (i-1)*w+1:i*w, 1,n);
        sely = clamp( (j-1)*w+1:j*w, 1,n);
        X = clamp(x + X0 + dX,1,n);
        Y = clamp(y + Y0 + dY,1,n);
        % extract patches
        P2 = M2(selx,sely);
        P1 = interp2( 1:n,1:n, M1, Y,X );
        % Compute best match and report its value.
        d = sum(sum( (P1-repmat(P2,[1 1 size(P1,3) size(P1,4)])).^2 ) );
        [tmp,I] = compute_min(d(:));
        F(selx,sely,1) = dx(I);
        F(selx,sely,2) = dy(I);
    end
end
