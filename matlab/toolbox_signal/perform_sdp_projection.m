function X = perform_sdp_projection(X)

% perform_sdp_projection - project on the cone of SDP matrices.
%
%   X = perform_sdp_projection(X);
%
%   Copyright (c) 2013 Gabriel Peyre

[V,D] = eig(X);
X = V*max(real(D),0)*V';

end