function [L,gA,gb,gx,X] = EvalNN(x,y, A,b,loss,rho)

% EvalNN - evaluate and compute the gradient of a fully connected network
%
%   [L,gA,gb,gx] = EvalNN(x,y, A,b,loss,rho);
%
%   Copyright (c) 2017 Gabriel Peyre

q = size(x,2); % #points
K = length(A);

%%% FORWARD %%%
X = {};
X{1} = x;
for k=1:K
    X{k+1} = rho.F( A{k}*X{k}+repmat(b{k},[1 q]) );
end

if isempty(y)
    L = X; 
    return;
end

L = loss.F(X{K+1},y);

%%% BACKWARD %%%
gx = loss.G(X{K+1},y); % nabla_X{k+1}
gA = {}; gb = {};
for k=K:-1:1
    M = rho.G(A{k}*X{k}+repmat(b{k},[1 q])) .* gx;
    % nabla_X{k} = [Df/dXk]^*( nabla_X{k+1} )
    gx = A{k}' * M;
    % nabla_A{k} = [Df/dAk]^*( nabla_X{k+1} )
    gA{k} =  M * X{k}';
    % nabla_b{k} = [Df/dbk]^*( nabla_X{k+1} )
    gb{k} =  sum(M,2);
end


end