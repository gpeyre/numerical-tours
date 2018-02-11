function [e,v,g] = compute_butterfly_neighbors(k, nj)

%   compute_butterfly_neighbors - compute local neighbors of a vertex
%
%   [e,v,g] = compute_butterfly_neighbors(k, nj);
%
%   This is for internal use.
%
% e are the 2 direct edge neighbors
% v are the 2 indirect neighbors
% g are the fare neighbors
%
%   You need to provide:
%       for e: vring, e2f 
%       for v: fring
%       for g: facej
%
%   Copyright (c) 2007 Gabriel Peyre

global vring e2f fring facej;

% find the 2 edges in the fine subdivition
vr = vring{k};
I = find(vr<=nj); 
e = vr(I);
% find the coarse faces associated to the edge e
f = [e2f(e(1),e(2)) e2f(e(2),e(1))];
% symmetrize for boundary faces ...
f(f==-1) = f(3 - find(f==-1));

if nargout>1
    F1 = mysetdiff(fring{f(1)}, f(2));
    F2 = mysetdiff(fring{f(2)}, f(1));
    % symmetrize for boundary faces
    F1 = [F1 repmat(f(1), 1, 3-length(F1))];
    F2 = [F2 repmat(f(2), 1, 3-length(F2))];
    v = [ mysetdiff( facej(:,f(1)), e )   mysetdiff( facej(:,f(2)), e ) ];
    if nargout>2
        d = [v, e];
        g = [setdiff( facej(:,F1(1)),d ), ...
            mysetdiff( facej(:,F1(2)),d ), ...
            mysetdiff( facej(:,F2(1)),d ), ...
            mysetdiff( facej(:,F2(2)),d ) ];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = mysetdiff(a,b)
% removed in a entries equal to entries in b
for s=b
    a(a==s) = [];
end