sel = randperm(size(P,1)); sel = sel(1:npts);
clf;
plot(P(sel,1), P(sel,2), '.');
axis([-1 1 -1 1]*5);
set_graphic_sizes([], 20);
title('Transformed domain');
