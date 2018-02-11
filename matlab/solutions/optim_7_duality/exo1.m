niter = 400;
Elist =[]; % primal
Flist = []; % dual
for i=1:niter
    u = proxG( u - gamma * nablaF(u), gamma );
    x = y - As(u);
    % 
    Elist(end+1) = E(x);
    Flist(end+1) = F(u);
    
end
Eplot = @(e)plot( log10(e(1:end/2)-min(e)) );
clf;
h = plot([Elist(7:end); -Flist(7:end)]'); axis tight;
set(h, 'LineWidth', 2);
legend('E', '-F');
