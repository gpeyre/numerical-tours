Energy = @(x)sqrt(sum(x(I).^2));
SoftAtten = @(x,gamma)max(0, 1-gamma./abs(x));
EnergyAtten = @(x,gamma)repmat(SoftAtten(Energy(x),gamma), [w*w 1] );
Flatten = @(x)x(:);
ProxG = @(x,gamma)accumarray(I(:), Flatten(EnergyAtten(x,gamma)), [N 1], @prod) .* x;
rProxG = @(x,gamma)2*ProxG(x,gamma)-x;
