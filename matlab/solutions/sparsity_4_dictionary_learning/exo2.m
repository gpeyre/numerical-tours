niter_dico = 100;
E = [];
tau = 1/norm(X*X');
for i=1:niter_dico
    R = D*X - Y;
    E(end+1) = sum(R(:).^2);
    D = ProjC( D + tau * (Y-D*X)*X' );
end
clf;
plot(log10(E(1:end/2)-min(E)));
axis tight;
