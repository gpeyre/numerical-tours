niter = 5;
u = ones(n,n,n);
err = [];
for i=1:niter
    Xrw = solvel1( A*diag(u(:)) ) .* u;
    u = (abs(Xrw).^(1-alpha)+epsilon);
    % 
    [tmp,xrw] = min( abs(Xrw-1), [], 3 );
    Xrw1 = encode(xrw); 
    err(end+1) = sum(A*Xrw1(:)~=1);
end
[tmp,xrw] = min( abs(Xrw-1), [], 3 );
Xrw1 = encode(xrw); 
disp(['Number of violated constraints: ' num2str(sum(A*Xrw1(:)~=1)) '.']);
