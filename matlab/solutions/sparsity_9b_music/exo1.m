% SVD
[U,S,V] = svd(MusicHankel(y),0); S = diag(S);
Ubot = U(:,N+1:end);
% d_y
d = Ubot'*exp(-2i*pi*(0:L-1)'*z(:)');
d = sum(abs(d).^2) / L;
% disp
clf; hold on;
plot(z, d, 'b');
stem(x0, 1+x0*0, 'r.');
axis([0 1 0 1]); box on;
