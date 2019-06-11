tau = .02; % step size
%
niter = 150*4;
q = 6;
ndisp = round( 1 + (niter/4-1)*linspace(0,1,q).^2 );
kdisp = 1;
%
f = zeros(n,1);
clf; 
E = [];
for it=1:niter
    % compute Laguerre cells and c-transform
    D = sum(Y.^2)' + sum(X.^2) - 2*Y'*X - f(:)';
    [fC,I] = min(D,[], 2);
    I = reshape(I, [p p]);
    E(it) = sum(f.*a) + sum(fC.*b(:));    
    % display
    if kdisp<=q && ndisp(kdisp)==it
        subplot(2,3,kdisp);
        hold on;
        imagesc(t,t,I); axis image; axis off;
        contour(t,t,I, -.5:n+.5, 'k', 'LineWidth', 2);
        colormap(Col);
        scatter( X(2,:), X(1,:), (1+(f-min(f))*10)*100, Col*.8, 'filled' );
        drawnow;
        kdisp = kdisp+1;
    end
    % sum area captured
    r = sum(sum( ( I==reshape(1:n,[1 1 n]) ) .* b, 1 ),2); r = r(:);
    % gradient ascent
    f = f + tau*( a-r(:) );
end
