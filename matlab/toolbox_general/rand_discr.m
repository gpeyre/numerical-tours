function sample = rand_discr(p, m)

% rand_discr - discrete random generator 
% 
%   y = rand_discr(p, n);
%
%   y is a random vector of length n drawn from 
%   a variable X such that
%       p(i) = Prob( X=i )
%
%   Copyright (c) 2004 Gabriel Peyré


if nargin<2
    m = 1;
end

% makes sure it sums to 1
p = p(:)' / sum(p(:));

n = length(p);

coin=rand(1,m);
cumprob=[0 cumsum(p)];
sample=zeros(1,m);
for j=1:n
  ind=find((coin>cumprob(j)) & (coin<=cumprob(j+1)));
  sample(ind)=j;
end
