% position
[x1,I] = sort(mod(angle(R),2*pi)/(2*pi));
R = R(I);
% Keep only the best N ones.
[~,I] = sort(abs(abs(R)-1));
x1 = x1(I(1:N));
R = R(I(1:N));
% Compute amplitude by solving a least square.
a1 = real(Phi(x1)\y);
% Display the recovered measure.
clf; hold on;
plot(z, f);
stem(x0,a0, 'r.');
stem(x1,a1, 'k--');
