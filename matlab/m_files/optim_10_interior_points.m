%%
% Interior points on the Lasso.

addpath('NewtonIP/');

n = 40; % #samples
p = 60; % dimension of the model
dotp = @(x,y)sum(x(:).*y(:));

I = round([.1 .25 .6 .9]*p);
w0 = zeros(p,1); w0(I) = [.8 -.6 .7 -.9];

X = randn(n,p);

y = X*w0 + .03*randn(n,1)*max(abs(X*w0)); % y=A*x0+noise
% y = randn(n,1);

lambda = max(abs(X'*y))/10; % regularization parameter

%% 
% The Lasso reads
%    min_w 1/2*|A*w-y|^2 + lambda*|w|_1
% which we re-write in primal as
%   min_{x=(u,v)>=0} 1/2*|A*(u-v)-y|^2 + lambda*<u+v,1>

h = @(w)1/2*norm(X*w-y)^2 + lambda*norm(w,1);

%%
% We consider the following general type of quadratic minimization under
% linear constraint
%       min_{x in R^d, C*x<=b} 1/2*<Q*x,x> + <x,s>

d = 2*p; m = 2*p;
A = -eye(m); b = zeros(m,1); 
Q = [X'*X, -X'*X; -X'*X, X'*X];
s = [lambda - X'*y; lambda + X'*y];

f = @(x)1/2*dotp(Q*x,x) + dotp(x,s);

% should be 0
w = randn(p,1);
x = [max(w,0);max(-w,0)];
disp( ['Should be 0: ', num2str(h(w) - (f(x)+ + 1/2*norm(y)^2))] );

            
%% 
% Newton using backtracking

x = rand(m,1); % initial point
t = 10; % not accurate 
t = 100; % medium accurate 
epsilon = 1e-10; % tolerance of newton
t = 1000; % very accurate 
[x,Flist,nbr_linesearch] = NewtonIP(Q,s,A,b,t,x, epsilon);

clf;
subplot(2,1,1);
plot(Flist);
subplot(2,1,2);
plot(nbr_linesearch);


clf; hold on;
stem(w0, 'k--');
stem(x(1:end/2), 'r', 'filled');
stem(-x(end/2+1:end), 'b', 'filled');
box on;

%%
% The dual of the lasso reads
%   min_{|X'*z|<=1} -lambda/2*|z|^2 + <z,y>
% with primal-dual relation that 
%   z = y-A*w  and   {i : |A'*z|=1}={i : w_i \neq 0}

A = [X';-X'];
b = ones(2*p,1);
Q = lambda/2*eye(n);
s = -y;

% run newton on the dual
z = zeros(n,1); % initial point
[z,Flist,nbr_linesearch] = NewtonIP(Q,s,A,b,t,z, epsilon);

w = x(1:end/2)-x(end/2+1:end);
eta = X'*z;
clf; hold on;
stem(.8*w/max(w), 'r', 'filled');
stem(eta, 'b', 'filled');
plot([eta*0+1, eta*0-1], 'k--');
axis([1 p -1.05 1.05]);
box on;
legend('X^T z', 'w');

%%
% Perform the IP method by progressively increasing the value of t.


mu = 20; % amplification factor
niter_ip_max = 100; % #iteration of IP

t = 1; % initial value
z = zeros(n,1); % initial point
step = 0;
flistIP = [];
for it=1:niter_ip_max
    [z,Flist,nbr_linesearch,flist] = NewtonIP(Q,s,A,b,t,z, epsilon);
    flistIP = [flistIP; flist(:)];
    t = t*mu;
    step(end+1) = step(end) + length(flist);  
    if length(b)/t < epsilon
        break
    end
end
  

clf; hold on;
Err = flistIP-min(flistIP);
plot(log10( Err/Err(1) ), 'LineWidth', 2);
for it=1:length(step)
    plot([step(it) step(it)], [-15 0], 'k--');
end
axis([0 step(end) -15 0]); box on;

%%
% test for different value of mu

mulist = [2 10 50 100 500 1000];
flistIP = {};
for k=1:length(mulist)
    mu = mulist(k);
    %
    t = 1; % initial value
    z = zeros(n,1); % initial point
    flistIP{k} = [];
    for it=1:niter_ip_max
        [z,Flist,nbr_linesearch,flist] = NewtonIP(Q,s,A,b,t,z, epsilon);
        flistIP{k} = [flistIP{k}; flist(:)];
        t = t*mu;
        if length(b)/t < epsilon
            break
        end
    end    
end

fmin = +Inf;
for k=1:length(mulist)
    fmin = min(fmin, min(flistIP{k}));
end
clf; hold on;
for k=1:length(mulist)
    s = (k-1)/(length(mulist)-1);
	Err = flistIP{k}-fmin;
    plot(log10( Err/Err(1) ), 'LineWidth', 2, 'color', [s 0 1-s]);
end
axis([0 step(end) -10 0]); box on;
% display if one use a single newton method
[~,~,~,flist] = NewtonIP(Q,s,A,b,1e5, zeros(n,1), epsilon);
plot(log10( (flist-fmin)/Err(1) ), 'k--', 'LineWidth', 2);


