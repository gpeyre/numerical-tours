function x = perform_saturation(x,tau,use_mad)

% perform_saturation - saturate a vector for better contrast
%
%   x = perform_saturation(x,tau,use_mad);
%
%   tau (around 1) is the saturation factor
%
%   copyright (c) 2007 Gabriel Peyre

if nargin<2
    tau = 1;
end
if nargin<3
    use_mad = 1;
end

tau = 2;
x = x-mean(x(:)); 
x = clamp( x/(2*mad(x(:))), -tau,tau);