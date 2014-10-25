r = .95*n/2;
p = 128; # number of points on the curve
theta = transpose( linspace(0, 2*pi, p + 1) )
theta = theta[0:-1]
gamma0 = n/2*(1+1j) +  r*(cos(theta) + 1j*sin(theta))
gamma = gamma0
clf;
imageplot(transpose(f))
cplot(gamma, 'r', 2)

