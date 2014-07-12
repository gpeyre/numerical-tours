switch setting
    case 1
        lambda_max = .005;
    case 2
        lambda_max = .01;
    otherwise
        error('Unknown setting 4');
end
lambda_list_TI = linspace(lambda_max/4,lambda_max,niter);
a = PsiS(y)*0;
errTI = [];
% Warmup
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI);
    a = SoftThresh( a + tau*PsiS(Phi(d)), lambda_list_TI(1)*tau );
end
% descent
for i=1:niter
    fTI = Psi(a);
    d = y-Phi(fTI);
    % step 
    a = SoftThresh( a + tau*PsiS(Phi(d)), lambda_list_TI(i)*tau );
    errTI(end+1) = snr(f0,fTI);
	if i>1 && errTI(i)>max(errTI(1:i-1))
        fBestTI = fTI;
    end
end
clf;
plot(lambda_list_TI,errTI); axis('tight');
