niter = 1000;
clf;
f1 = f;
disp_list = [1 3 10 100]; q = 1;
for i=1:niter    
    [Theta,~] = qr(randn(d));    
    f1 = (1-tau)*f1 + tau * Theta * P(Theta'*f1, Theta'*g);
    if q<=4 && i==disp_list(q)
        t = (q-1)/3;
        subplot(2,2,q);
        hold on;
        plotp(f1, [t 0 1-t]);
        axis('off'); axis('equal');
        q = q+1;
    end
end
