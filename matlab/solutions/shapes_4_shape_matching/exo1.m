[Y,X] = meshgrid(1:n,1:n);
d = [];
for k=1:nrad
    r = rlist(k);    
    if 0
        Mh = perform_blurring(M1,r);
    else
        x = -ceil(r):ceil(r);
        [b,a] = meshgrid(x,x);
        h = double( a.^2+b.^2<=r^2 );
        h = h/sum(h(:));
        Mh = perform_convolution(M,h);
    end
    [Y,X] = meshgrid(1:n,1:n);
    for i=1:2
        D{i}(k,:) = interp2(Y,X,Mh{i},bound{i}(2,:), bound{i}(1,:));
        % I = round(bound{i}(1,:)) + round(bound{i}(2,:)-1)*n;
        % D{i}(k,:) = Mh{i}(I);
    end
end
clf;
sel = 1:round(nbound/8);
plot(sel, D{1}([1 nrad/2 nrad],sel)', '-');
title('Some features (zoom)');
axis tight;
