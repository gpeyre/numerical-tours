function y = keep_biggest(x,n)

% keep_biggest - keep only the n biggest coef, 
%   set the rest to zero.
%
% y = keep_biggest(x,n);
%
%   Copyright (c) 2004 Gabriel Peyré

if iscell(x)
    % concat cells array
    xx = [];
    for i=1:length(x)
        xx = [xx;x{i}(:)];
    end
    yy = keep_biggest(xx,n);
    % recreate cells
    for i=1:length(x)
        y{i} = yy( 1:prod( size(x{i}) ) );
        y{i} = reshape(y{i}, size(x{i}));
        yy( 1:prod( size(x{i}) ) ) = [];
    end
    return;
end

n = round(n);
y = x;
[tmp,I] = sort(abs(y(:)));
y( I(1:end-n) ) = 0;