niter = 60;
q = 6;
ndisp = round( 1 + (niter/4-1)*linspace(0,1,q).^2 );
kdisp = 1;
%
X1 = X;
clf; 
E = [];
for it=1:niter
    % compute Laguerre cells and c-transform
    D = sum(Y.^2)' + sum(X1.^2) - 2*Y'*X1;
    [fC,I] = min(D,[], 2);
    I = reshape(I, [p p]);
    E(it) = sum(fC.*b(:));    
    % display
    if kdisp<=q && ndisp(kdisp)==it
        subplot(2,3,kdisp);
        hold on;
        imagesc(t,t,I); axis image; axis off;
        contour(t,t,I, -.5:n+.5, 'k', 'LineWidth', 2);
        colormap(Col);
        scatter( X1(2,:), X1(1,:), 100, 'k', 'filled' );
        drawnow;
        kdisp = kdisp+1;
    end
    % update barycenter    
    A = ( I==reshape(1:n,[1 1 n]) ) .* b;
    B = ( I==reshape(1:n,[1 1 n]) ) .* b .* ( U+1i*V );
    X1 = sum(sum(B,1),2) ./ sum(sum(A,1),2); 
    X1 = [real(X1(:))';imag(X1(:))'];
end
