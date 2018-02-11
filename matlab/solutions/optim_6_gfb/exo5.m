nIter = 500;
lambda = bestLambda;
mu = 1e-3;
x = zeros(N,N); % initialization minimizer
z = zeros(N,N,n); % initialization auxiliary variables
for it=1:nIter
	forward = 2*x - gamma*nablaF(x);
	for i=1:4
		z(:,:,i) = z(:,:,i) + proxGi(forward-z(:,:,i),n*gamma*lambda,i) - x;
	end
	z(:,:,5) = z(:,:,5) + proxG5(forward-z(:,:,5),n*gamma*mu) - x;
	x = mean( z, 3 );
end
clf
imageplot(x)
title( sprintf( '\\lambda=%.1e; \\mu=%.1e; SNR=%.2fdB', lambda, mu, snr(f0,x) ) );
