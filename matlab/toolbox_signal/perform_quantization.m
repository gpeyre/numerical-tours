function [y, y_quant] = perform_quantization(x, T, dir)

% perform_quantization - perform a quantization of the signal with double zero bin.
%
% If x is a double array (or a cell collection of arrays) 
% to be quantized, use
%   [y, y_quant] = perform_quantization(x, T, 1);
% y is the quantified value (real value)
% y_quant is the signed token representing each entry of y
%   (this is an integer).
% If y_quant is an int32
%   x = perform_quantization(y_quant, T, -1);
%
%   The quantizer is defined by y=Q_T(x) where:
%       Q_T(x) = 0    if  |x|<T
%       Q_T(x) = sign(x) * ([x/T]+0.5)*T      where [.]=floor
%   (i.e. a nearly uniform quantizer with twice larger zero bin).
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(x)
    if nargin<3
        if isinteger(x{1})
            dir = -1;
        else
            dir = 1;
        end
    end
    % recurse on each entry of the cell array
    y = {}; y_quant = {};
    for i=1:length(x)
        [y{i}, y_quant{i}] = perform_quantization(x{i}, T, dir);
    end
    return;
end


if nargin<3
    dir = 1;
end
    
if dir==-1
    % perform un-transformation
    x = double(x);
    y = sign(x) .* (abs(x)+0.5) * T;
    return;
end

% perform quantization
I = find(abs(x)<T);
q = floor(abs(x)/T);

y_quant = sign(x).*q;

y = sign(x).*(q+0.5)*T;
y(I) = 0;