d = D(:,:,n-delta);
[tmp,I] = min(d(:));
[x,y] = ind2sub([n n],I(1));
end_point = [x;y;n-delta];
