niter = 1000;
clf;
f1 = f;
q = 1; disp_list = [3 10 100 niter];
for i=1:niter    
    [Theta,~] = qr(randn(d));    
    f1 = (1-tau)*f1 + tau * Theta * P(Theta'*f1, Theta'*g);
    if q<=4 && i==disp_list(q)
        subplot(2,2,q);
        F1 = reshape(f1', [n n 3]);
        imageplot(F1);
        q = q+1;
    end
end
%
F1 = reshape(f1', [n n 3]);
