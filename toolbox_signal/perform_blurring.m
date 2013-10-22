function M = perform_blurring(M, sigma, options)

% perform_blurring - gaussian blurs an image
%
%   M = perform_blurring(M, sigma, options);
%
%   M is the original data
%   sigma is the width of blurs (in pixels)
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;

if iscell(M)
    for i=1:length(M)
        M{i} = perform_blurring(M{i}, sigma, options);
    end
    return;
end

if sigma==0
    return;
end

if size(M,3)>1
    for i=1:size(M,3)
        for j=1:size(M,4)
            M(:,:,i,j) = perform_blurring(M(:,:,i,j), sigma, options);
        end
    end
    return;
end

n = max(size(M));

eta = 4;
p = round((sigma*eta)/2)*2+1;
p = min(p,round(n/2)*2-1);

A = [1 1];
if size(M,1)==1 || size(M,2)==1
    A = 1; % 1D.
end
h = compute_gaussian_filter(p*A,sigma/(4*n),n*A);
M = perform_convolution(M, h, options);