glist = [.1 .01 .005 .001 ];
niter = 300;
clf;
for k=1:length(glist)
    gamma = glist(k);
    xi = exp(-C/gamma);
    b = ones(N(2),1); 
    for i=1:niter
        a = p ./ (xi*b);
        b = q ./ (xi'*a); 
    end
    Pi = diag(a)*xi*diag(b);
    imageplot( clamp(Pi,0,min(1./N)*.3) , ['\gamma=' num2str(gamma)], 2,2,k);
end
