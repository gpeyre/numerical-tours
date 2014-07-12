r = .95*n/2;
p = 128; % number of points on the curve
theta = linspace(0,2*pi,p+1)'; theta(end) = [];
gamma0 = n/2*(1+1i) +  r*(cos(theta) + 1i*sin(theta));
gamma = gamma0;
clf; hold on;
imageplot(f);
h = plot(imag(gamma([1:end 1])),real(gamma([1:end 1])), 'r');
set(h, 'LineWidth', 2);
