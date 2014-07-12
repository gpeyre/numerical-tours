mlist = N*[10 5 2 1];
clf;
for i=1:length(mlist)
    subplot(4,1,i);
    m = mlist(i);
    S1 = sparse(Gamma(S,tau(m)));
    plot( PsiS( S1*Psi(f) )  );
    axis tight;
    title(['m/N=' num2str(m/N)]);
end
