niter = 5000;
C = zeros(n,n,k);
for i=1:k
    u = zeros(k+1,1);
    u(i) = 1;
    if i==1
        u(end)=1;
    end
    u = interp1(0:k,u, abscur);
    Ci = zeros(n);
    for it=1:niter
        Ci = ( Ci(sel1,:) + Ci(:,sel1) + Ci(sel2,:) + Ci(:,sel2) )/4;
        Ci(I) = u;
    end
    C(:,:,i) = Ci;
end
