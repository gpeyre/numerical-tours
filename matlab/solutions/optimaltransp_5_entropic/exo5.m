m = 5; 
[T,S] = meshgrid(linspace(0,1,m), linspace(0,1,m));
clf;
niter = 1000;
for j=1:m^2
    % weights
    lambda = [S(j)*T(j), (1-S(j))*T(j), S(j)*(1-T(j)), (1-S(j))*(1-T(j))];
    % computation
    b = ones(N,N,K); a = b;
    for i=1:niter
        for k=1:K
            a(:,:,k) = P(:,:,k) ./ xi(b(:,:,k));
        end
        q = zeros(N);
        for k=1:K
            q = q + lambda(k) * log( max(1e-19, b(:,:,k) .* xi(a(:,:,k)) ) );
        end
        q = exp(q);
        for k=1:K
            b(:,:,k) = q ./ xi(a(:,:,k));
        end
    end
    % display
    subplot(m,m,j);
    imageplot(q); drawnow;
end
