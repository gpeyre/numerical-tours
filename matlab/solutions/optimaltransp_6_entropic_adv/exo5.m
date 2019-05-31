A = eye(2); h = zeros(2,1);
tau_A = 1; % step size for the descent 
tau_h = .15;
giter = 40; % iter for the gradient descent
ndisp = round(linspace(1,giter,6)); kdisp = 1;% ndisp(1) = [];
f = zeros(n,1);
clf;
for j=1:giter
    x = A*z+h;
    C = distmat(x,y); 
    for it=1:niter
        g = mina(C-f,epsilon);
        f = minb(C-g,epsilon);
    end
    P = a .* exp((f+g-C)/epsilon) .* b;
    % gradient with respect to positions
    v = a' .* x - ( y*P' );
    % gradient with respect to parameters
    nabla_A = v*z';
    nabla_h = sum(v,2);
    A = A - tau_A * nabla_A;
    h = h - tau_h * nabla_h;    
    % drawing
    if ndisp(kdisp)==j
        subplot(2,3,kdisp);
        s = (j-1)/giter;
        hold on;
        plotp(y, 'r');
        plotp(x, 'b');
        axis equal; axis tight; axis off;
        drawnow;
        kdisp = kdisp+1;
    end
end
%EXO