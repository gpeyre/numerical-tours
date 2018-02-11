function [h,x] = compute_histogram(M, options)

% compute_histogram - compute the (symmetric) histogram of a vector after threshold.
%
% [h,x] = compute_histogram(M,options);
%
%   If M contains only entry in 1...p, where p=max(M), then h is of size p
%       h(k) = #{i \ M(i)=k} / length(M)
%   If M contains only entry in -p...p, where p=max(abs(M)), then h is of size 2*p+1
%       h(k) = #{i \ M(i)=k-p-1} / length(M)
%   (you can also force this by setting options.symmetrize=1).
%
%   For real valued signals, you can add thresholding by setting options.T
%
%   Copyright (c) 2006 Gabriel Peyré

if nargin>1 && ~isstruct(options)
    error('Valid syntax is "h = compute_histogram(M,options);".');
end

options.null = 0;
if isfield(options, 'avoid_zeros')
    avoid_zeros = options.avoid_zeros;
else
    avoid_zeros = 1;
end
if isfield(options, 'T')
    T = options.T;
else
    T = -1;
end

if length(T)>1
	% build a cell array of histograms
	h = {};
	for i=1:length(T)
		h{end+1} = compute_histogram(M,T(i));
	end
	return;
end

M = M(:);

if isempty(M)
    h = []; x = [];
    return;
end

% quantize
if T>0
	[Mw, M] = perform_quantization(M, T);
    M = double(M);
end

if isfield(options, 'symmetrize')
    symmetrize = options.symmetrize;
    if ~symmetrize && ~isempty(find(M<=0))
        error('Negative entry do not work with options.symmetrize=0');
    end
else
    symmetrize = ~isempty(find(M<=0));
end    

% bounds of histogram
hp = max(abs(M(:)));
hpmax = 10000;
if hp>hpmax
    error(['Number of entries in histogram is > ' num2str(hpmax) '.']);
end
hp = min(double(hp),hpmax);   % avoid too big histogramms

% compute histrogram
if symmetrize
    h = zeros(2*hp+1,1);
    x = -hp:hp;
    for p = x
        h(p+hp+1) = sum(M==p);
    end
else
    h = zeros(hp,1);
    x = 1:hp;
    for p = x
        h(p) = sum(M==p);
    end 
end

% avoid zero value
if avoid_zeros
    I = find(h==0);
    J = find(h>0);
    h(I) = mmin(h(J))/1e12;
end

% normalize to sum to 1
h = h/sum(h);

function y = mmin(x)
y = min(x(:));