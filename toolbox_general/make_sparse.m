function A = make_sparse(i,j,x, n, m)

% make_sparse - synonymous with sparse(i,j,x)


if nargin<3
    A = sparse(i,i);
    return;
end

if nargin<4
    n = max(i);
end
if nargin<5
    m = max(j);
end

A = sparse( i, j, x, n, m );