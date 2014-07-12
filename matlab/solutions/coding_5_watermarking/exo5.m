niter = 4000;
c = [];
for i=1:niter
    w = randn(P,1);
    rho = sqrt(P)/norm(w.*abs(x0)) * 10^( -psnr_embedding/20 );
    x = x0 + rho*abs(x0).*w;  
    c(end+1) = C(A(x),w);
end
clf;
hist(c, 30);
