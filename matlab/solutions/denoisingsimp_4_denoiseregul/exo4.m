niter = 800;
lambda_list = linspace(0,.25,niter);
tau = 2 / ( 1 + max(lambda) * 8 / epsilon);
fTV = y;
energy = [];
for i=1:niter
    lambda = lambda_list(i);    
    Gr = grad(fTV);
    d = sqrt(sum3(Gr.^2,3));
    G0 = -div( Gr ./ repmat( sqrt( epsilon^2 + d.^2 ) , [1 1 2]) );
    G = fTV-y+lambda*G0;
    deps = sqrt( epsilon^2 + d.^2 );
    fTV = fTV - tau*G;
    err(i) = snr(f0,fTV);
    if i>1
        if err(i) > max(err(1:i-1))
            fTV0 = fTV;
        end
    end
end
clf;
plot(lambda_list, err); axis('tight');
set_label('\lambda', 'SNR');
