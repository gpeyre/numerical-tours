function f1 = perform_curve_subdivision(f, nsub, options)

% perform_curve_subdivision - perform subdivision
%
%   f1 = perform_curve_subdivision(f, nsub, options);
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if nsub<=0
    f1 = f; return;
end
if size(f,1)>size(f,2)
    f = f';
end
if size(f,1)>1
    for k=1:size(f,1)
        f1(k,:) = perform_curve_subdivision(f(k,:), nsub, options);
    end
    return;
end

h = getoptions(options, 'h', [1 4 6 4 1]);
h = 2*h/sum(h(:));

if not(isfield(options, 'bound'))
    options.bound = 'per';
end

f = perform_curve_subdivision(f, nsub-1, options);

f = f(:)';
n = length(f);
f1 = zeros(1,n*2);
f1(1:2:end) = f;
f1 = perform_convolution(f1, h, options);
f1 = f1(:)';