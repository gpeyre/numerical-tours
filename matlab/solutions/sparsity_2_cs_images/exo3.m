ProxF = @(x,gamma)x + PhiS(y-Phi(x)); 
ProxG = @(x,gamma)max(0,1-gamma./max(1e-15,abs(x))).*x;
rProxF = @(x,gamma)2*ProxF(x,gamma)-x;
rProxG = @(x,gamma)2*ProxG(x,gamma)-x;
