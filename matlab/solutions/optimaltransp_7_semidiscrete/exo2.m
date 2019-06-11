niter = 2000;
q = 6;
ndisp = round( 1 + (niter/2-1)*linspace(0,1,q).^2 );
kdisp = 1;
%
f = zeros(n,1);
clf; 
E = [];
for it=1:niter
    % sample
    i = (rand<W(1))+1; % select one of the two Gaussian
    y = [S(i) * randn + Z{i}(1);S(i) * randn + Z{i}(2)];
    % detect Laguerre cell where y is
    [~,i] = min( sum(y.^2)' + sum(X.^2) - 2*y'*X - f(:)' );
    % gradient
    nablaEy = a; nablaEy(i) = nablaEy(i) - 1; 
    % gradient ascent
    l0 = 10; % warmup phase.
    tau = .1/( 1 + it/l0 );
    f = f + tau*nablaEy;
    
    if 1 % mod(it,10)==1
        % compute Laguerre cells and c-transform
        D = sum(Y.^2)' + sum(X.^2) - 2*Y'*X - f(:)';
        [fC,I] = min(D,[], 2);
        I = reshape(I, [p p]);
        E(end+1) = sum(f.*a) + sum(fC.*b(:));
    end
    
    % display
    if kdisp<=q && ndisp(kdisp)==it
        D = sum(Y.^2)' + sum(X.^2) - 2*Y'*X - f(:)';
        [fC,I] = min(D,[], 2);
        I = reshape(I, [p p]);
        %
        subplot(2,3,kdisp);
        hold on;
        imagesc(t,t,I); axis image; axis off;
        contour(t,t,I, -.5:n+.5, 'k', 'LineWidth', 2);
        colormap(Col);
        scatter( X(2,:), X(1,:), (1+(f-min(f))*10)*100, Col*.8, 'filled' );
        drawnow;
        kdisp = kdisp+1;
    end
end
