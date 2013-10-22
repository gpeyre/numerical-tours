function [v,w] = perform_haar_graph(v, A, dir, options)

% perform_haar_graph - perform the haar transform on graph
%
% Forward transform:
%   [vwav,w] = perform_haar_graph(v, A, +1, options);
% Backward transform:
%   [v,w] = perform_haar_graph(vwav, A, -1, options);
%
%   A is the (sparse) adjacency matrix of the graph.
%   v is the signal.
%   vwav is its wavelet transform coefficients (orthogonal decomposition).
%   w is the weight (size of the Haar function)
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
verb = getoptions(options, 'verb', 1);

v = v(:);
n = length(v);
A = sparse(A);
n1 = n-2;

%% First simulate the grouping process %%
ilist = zeros(n1,1);
jlist = zeros(n1,1);
wilist = zeros(n1,1);
wjlist = zeros(n1,1);
w = ones(n,1);
% find (directed) edges
[I,J] = find(A);
a = find(I<J); I = I(a); J = J(a);
active = ones(n,1);
for iter=1:n1
    if isempty(I) || isempty(J)
        break;
    end
    if verb
        progressbar(iter,n1);
    end
    % find position with minimum weight
    [tmp,k] = min( w(I)+w(J) );
    i = I(k); j = J(k);
    ilist(iter) = i; jlist(iter) = j;
    wilist(iter) = w(i); wjlist(iter) = w(j);
    % new weight
    w(i) = w(i)+w(j);
    % redirect connexions of node j on node i
    I(I==j) = i; J(J==j) = i;
    a = find(I<J); I = I(a); J = J(a); 
    active(j) = 0;
end
if verb
    nmax = iter-1;
    progressbar(n1,n1);
end
% remaining coarse coefficient
icoarse = find(active);
%[icoarse,J] = I;
wcoarse = sqrt(w(icoarse));

if dir==+1
    %% Forward transform %%
    for iter=1:nmax
        i = ilist(iter);
        j = jlist(iter);
        wi = wilist(iter);
        wj = wjlist(iter);
        % detail coefficient
        d = (v(i)-v(j)) * sqrt( wi*wj/(wi+wj) );
        % average
        v(i) = ( v(i)*wi + v(j)*wj )/(wi+wj);
        v(j) = d;
    end
    % normalization of coarse scale
    v(icoarse) = v(icoarse) .* wcoarse;
else
    %% Backward transform %%
    % un - normalization of coarse scale
    v(icoarse) = v(icoarse) ./ wcoarse;
    for iter=nmax:-1:1
        i = ilist(iter);
        j = jlist(iter);
        wi = wilist(iter);
        wj = wjlist(iter);
        % detail coefficient
        d = v(j)*sqrt( (wi+wj)/(wi*wj) );
        m = v(i)*(wi+wj);
        v(j) = (m-wi*d)/( wi+wj );
        v(i) = d+v(j);
    end
end

return;

%% OLD CODE %%

for iter=1:n1
    if verb
        progressbar(iter,n1);
    end
    % find edges
    [I,J] = find(A);
    % find position with minimum weight
    [tmp,k] = min( w(I)+w(J) );
    i = I(k); j = J(k);
    ilist(iter) = i; jlist(iter) = j;
    wilist(iter) = w(i); wjlist(iter) = w(j);
    % new weight
    w(i) = w(i)+w(j);
    % redirect connexion of node j on node i
    A(i,:) = A(i,:) | A(j,:);
    A(:,i) = A(:,i) | A(:,j);
    A(i,i) = 0;
    A(j,:) = 0; A(:,j) = 0;
end