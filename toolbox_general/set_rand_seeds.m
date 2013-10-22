function set_rand_seeds(a,b)

% set_rand_seeds - initialize rand and randn

if nargin<1
    a = 123456;
end
if nargin<1
    b = 123456;
end
randn('state', a);
rand('state', b);