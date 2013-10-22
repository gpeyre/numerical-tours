function h = plot_levelset(M,t,A)

% plot_levelset - display the level set of an image
%
%   plot_levelset(M,t,A);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<2
    t = 0;
end
if nargin<3
    A = M;
end

hold on;
imageplot(A);
[c,h] = contour(M, [t t], 'r');
set(h, 'LineWidth', 2);
axis ij;
hold off;