niter = 500;
displist = round([.1 .2 .5 1]*niter); idisp = 1;
clf;
Ci = zeros(n);
for it=1:niter
    Ci = ( Ci(sel1,:) + Ci(:,sel1) + Ci(sel2,:) + Ci(:,sel2) )/4;
    Ci(I) = u;
    if it==displist(idisp)
        c = Ci + 1e-3;
        c(S==1) = 0;
        B = display_shape_function(c');
        subplot(2,2,idisp);
        hold on;
        t = linspace(0,1,n);
        imagesc(t,t,B); axis('image'); axis('off');
        colormap jet(256);
        h = plot(V(1,[1:end 1]), V(2,[1:end 1]), 'r.-');
        set(h, 'LineWidth', lw); set(h, 'MarkerSize', ms);
        idisp = idisp+1;
    end
end
C(:,:,i) = Ci;
