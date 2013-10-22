function [v1,v2] = compute_levelset_mesh(vertex,face,f,tau,options)

% compute_levelset_mesh - compute level set curve
%
%   [v1,v2] = compute_levelset_mesh(vertex,face,f,tau,options);
%
%   v1(:,i),v2(:,i) is a segment of a levelt set of point f(x)=tau on the
%   mesh.
%
%   Copyrigh (c) 2007 Gabriel Peyre


if length(tau)>1
    v1 = [];
    v2 = [];
    for i=1:length(tau)
        [w1,w2] = compute_levelset_mesh(vertex,face,f,tau(i),options);
        v1 = [v1, w1];
        v2 = [v2, w2];
    end
    return;
end

f = f-tau;

s = sum( f(face)>0 );
I = find( s<3 & s>0 );
m = length(I);

% find intersection points
t = f(face(:,I))./( f(face(:,I))-f(face([2 3 1],I)) );
t(t<0 | t>1) = Inf;
[tmp,J1] = min(t,[],1); t1 = t(J1+(0:m-1)*3);
t(J1+(0:m-1)*3) = Inf;
[tmp,J2] = min(t,[],1); t2 = t(J2+(0:m-1)*3);
A1 = face(J1+(I-1)*3);
A2 = face(mod(J1,3)+1+(I-1)*3);

B1 = face(J2+(I-1)*3);
B2 = face(mod(J2,3)+1+(I-1)*3);

v1 = vertex(:,A1).*repmat(1-t1,[3 1]) + vertex(:,A2).*repmat(t1,[3 1]);
v2 = vertex(:,B1).*repmat(1-t2,[3 1]) + vertex(:,B2).*repmat(t2,[3 1]);