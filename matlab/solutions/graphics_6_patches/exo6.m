noffs = 1;
niter = 6;
% initialization
M = M0;
I = find(mask==1);
M(I) = rand(length(I),1);
clf;
for it=1:niter
    % randomize the offsets
    sel = randperm(w*w); sel = sel(1:noffs);
    OffX = dX(sel); OffY = dY(sel);
    % project
    M1 = zeros(n);
    for j=1:noffs
        ofx = OffX(j);
        ofy = OffY(j);
        Xs = mod(X+ofx-1, n)+1;
        Ys = mod(Y+ofy-1, n)+1;
        P = M(Xs + (Ys-1)*n);
        for i=1:p*p
            if sum(sum(Pmask(:,:,i)))>0
                d = sum(sum( (P0 - repmat(P(:,:,i), [1 1 q])).^2 ) );
                [tmp,s] = min(d);
                P(:,:,i) = P0(:,:,s);
            end
        end
        M1(Xs + (Ys-1)*n) = M1(Xs + (Ys-1)*n) + P;
    end
    M = M1 / noffs;
    M(mask==0) = M0(mask==0);
    imageplot(M,'', 2, ceil(niter/2), it);
end
