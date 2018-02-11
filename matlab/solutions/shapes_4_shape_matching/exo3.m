% metric
W = rescale(C, .001,1);
% distance
[D,S] = perform_fast_marching(1./W, [1;1]);
% path
gpath = compute_geodesic(D,size(D));
