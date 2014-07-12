compute_w = @(I,c)sum( (I - repmat(reshape(c,[1 1 3]), [n n 1])).^2, 3);
w0 = compute_w(I,c0);
w1 = compute_w(I,c1);
Omega = w0<w1;
display_segmentation = @(u)repmat(reshape(c0,[1 1 3]), [n n 1]) .* repmat(u>.5, [1 1 3]) + ...
        repmat(reshape(c1,[1 1 3]), [n n 1]) .* repmat(u<.5, [1 1 3]);
clf; 
imageplot(display_segmentation(Omega));
