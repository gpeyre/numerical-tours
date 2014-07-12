p = [[96;72], [230;175], [50;125], [95;180]];
clf;
for i=1:1
    [f1,S] = perform_fast_marching(1./W, p(:,i));
    f1 = -f1*n;
    subplot(1,1,i);
    hold on;
    surf(f1);
    h = plot3(p(2,i), p(1,i), f1(p(1,i),p(2,i)), 'r.');
    set(h, 'MarkerSize', 30);
    colormap(gray(256));
    shading interp;
    axis('equal');
    view(110,45);
    axis('off');
    camlight;
end
