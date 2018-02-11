delta = zeros(N,N);
sel = X(1,:) + (X(2,:)-1)*q;
for i=1:N
    U = geod(X(:,i));
    delta(:,i) = U(sel);
end
delta = (delta+delta')/2;
