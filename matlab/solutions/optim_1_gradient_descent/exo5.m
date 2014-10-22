epsilon = 200;
tau = 1.8/( 8/epsilon );
% tau = tau * 100;
E = [];
niter = 150;
ndisp = round(linspace(1,niter, 5)); ndisp(1) = [];
x = y;
clf; q = 1;
for i=1:niter
    E(i) = J(x,epsilon);
    if i>1 && E(i)>E(i-1)
       tau = tau*.8;
    end
    x = x - tau * GradJ(x,epsilon);
    x = ProjH(x,y);
    if i==ndisp(q)
        subplot(2,2,q);
        imageplot(x);
        q = q+1;
    end
end
