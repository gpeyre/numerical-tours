slist = [4 6 10 15];
clf;
for i=1:length(slist)
    sigma = slist(i);
    %
    g = grad( blur(f0,sigma) );
    h = hessian( blur(f0,sigma) );
    a = h(:,:,1:2).*repmat(g(:,:,1), [1 1 2]) + ...
        h(:,:,2:3).*repmat(g(:,:,2), [1 1 2]);
    %
    subplot( 2,2,i );
    plot_levelset( sum(a.*g, 3) ,0,f0);
    title(['\sigma=' num2str(sigma)]);
end
%EXO