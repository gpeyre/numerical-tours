z = x; % initialization
tau = .05; % step size for the descent 
giter = 20; % iter for the gradient descent
ndisp = round(linspace(1,giter,6)); kdisp = 1;% ndisp(1) = [];
f = zeros(n,1);
clf;
for j=1:giter
    % Sinkhorn
    C = distmat(z,y);
    for it=1:niter
        g = mina(C-f,epsilon);
        f = minb(C-g,epsilon);
    end
    P = a .* exp((f+g-C)/epsilon) .* b;
    % gradient
    G = z - ( y*P' ) ./ a';
    z = z - tau * G;
    % drawing
    if ndisp(kdisp)==j
        subplot(2,3,kdisp);
        s = (j-1)/giter;
        hold on;
        % plotp(x, 'b');
        % plotp(y, 'r');
        plotp(z, [s 0 1-s]);
        axis equal; axis tight; axis off;
        drawnow;
        kdisp = kdisp+1;
    end
end
