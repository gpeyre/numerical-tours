C = zeros(n,n,k);
D = zeros(n,n,k);
for i=1:k
    j = mod(i,k)+1;
    vi = V(:,i);
    vj = V(:,j);
    ni = N(:,i);
    %
    a = repmat(vj - vi,[1 n^2]);
    b = repmat(vi,[1 n^2]) - W;
    Q = sum(a.^2); s = sum(b.^2); R = 2*sum(a.*b);
    na = sqrt(sum(a.^2));
    BA = na .* sum( b .* repmat(ni,[1 n^2]) );
    SRT = sqrt( 4*s.*Q - R.^2 );
    L0 = log(s); L1 = log(s+Q+R);
    A0 = atan(      R ./SRT ) ./ SRT;
    A1 = atan( (2*Q+R)./SRT ) ./ SRT;
    A10 = A1 - A0;
    L10 = L1 - L0;
    %
    d = - na .* ( (4*s-(R.^2)./Q) .* A10 + R./(2*Q).*L10 + L1 - 2 )  / (4*pi) ;
    cj = + BA .* ( L10./(2*Q) - A10 .* (  R./Q) ) / (2*pi);
    ci = - BA .* ( L10./(2*Q) - A10 .* (2+R./Q) ) / (2*pi);
    D(:,:,i) = reshape(d,n,n);
    C(:,:,i) = C(:,:,i) + reshape(ci,n,n);
    C(:,:,j) = C(:,:,j) + reshape(cj,n,n);
end
