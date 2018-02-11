nIter = 100;
ObjList = zeros(1,nIter);
for it=1:nIter
	ObjList(it) = E(x,lambda);
	forward = 2*x - gamma*nablaF(x);
	for i=1:4
		z(:,:,i) = z(:,:,i) + proxGi(forward-z(:,:,i),n*gamma*lambda,i) - x;
	end
	x = mean( z, 3 );
end
clf
h = plot(ObjList);
% h = plot(log10(ObjList(1:round(end*.7))-min(ObjList)));
set(h, 'LineWidth', 2);
xlabel( 't' );
% title( 'log_{10}(E(x_t,\lambda) - E(x^{*},\lambda))' );
title( 'E(x_t,\lambda)' );
axis tight;
