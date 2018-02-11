if 0
qx = 10; qv = 10;
sxlist = linspace(.5,2*1.5,qx);
svlist = linspace(.05,.3,qv);
err = [];
for ix=1:qx
    for iv=1:qv
        sx = sxlist(ix); sv = svlist(iv);
        f1 = bilateral_lin(f,sx,sv);
        err(ix,iv) = snr(f0,f1);
        if err(ix,iv)==max(err(:))
            fOpt = f1;
        end
    end    
end
clf;
surf(sxlist,svlist,err');
axis('tight');
view(60,40);
colormap(jet(256));
xlabel('\sigma_x');
ylabel('\sigma_v');
else
    fOpt = f;
end
