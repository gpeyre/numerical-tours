nIter = 100;
lambdaListLen = length(lambdaList);
SNRlist = zeros(1,lambdaListLen);
SNRmax = 0;
bestLambda = 0;
for l=1:lambdaListLen
	lambda = lambdaList(l);
	x = zeros(N,N); % initialization minimizer
	z = zeros(N,N,n); % initialization auxiliary variables
	for it=1:nIter
        	forward = 2*x - gamma*nablaF(x);
		for i=1:4
		    z(:,:,i) = z(:,:,i) + proxGi(forward-z(:,:,i),n*gamma*lambda,i) - x;
		end
		x = mean( z, 3 );
	end
	SNRlist(l) = snr(f0,x);
	if SNRlist(l) > SNRmax
		recov = x;
		SNRmax = SNRlist(l);
		bestLambda = lambda;
	end
end
clf
h = semilogx( lambdaList, SNRlist );
set(h, 'LineWidth', 2);
xlabel( '\lambda' )
ylabel( 'SNR' )
axis tight
