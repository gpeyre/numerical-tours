niter = 800;
b = ones(N,N,K); a = b;
mynorm = @(x)norm(x);
Err_q = [];
for i=1:niter
    Err_q(i) = 0;
    for k=1:K
        Err_q(i) = Err_q(i) + mynorm( a(:,:,k) .* xi(b(:,:,k)) - P(:,:,k) )/mynorm(P(:,:,k));
        a(:,:,k) = P(:,:,k) ./ xi(b(:,:,k));
    end
    q = zeros(N);
    for k=1:K
        q = q + lambda(k) * log( max(1e-19, b(:,:,k) .* xi(a(:,:,k)) ) );
    end
    q = exp(q);
    % clf; imageplot(q); drawnow;
    for k=1:K
        b(:,:,k) = q ./ xi(a(:,:,k));
    end
end
clf;
plot(log(Err_q)); axis tight;
