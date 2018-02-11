function Q = perform_sos_projection(Q)

% perform_sos_projection - project on the linear constraint for SOS polynomials
%
%   Q = perform_sos_projection(Q);
%
%   The constraint is that each non-leading diagonal should sum to 0 and
%   the main diagonal should sum to 1.
%
%   Copyright (c) 2013 Gabriel Peyre

n = size(Q,1);
% index set
U = reshape(1:n^2, [n n]);
for i=-n+1:n+1
    I = diag(U,i);
    Q(I) = Q(I) - mean(Q(I));
end
% add back id/n so that diag sum to 1.
Q = Q + eye(n)/n;

end