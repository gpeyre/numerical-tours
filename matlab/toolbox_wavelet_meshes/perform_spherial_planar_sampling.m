function posw = perform_spherial_planar_sampling(pos_sphere, sampling_type, options)

% perform_spherial_planar_sampling - project sampling location from sphere to a square
%
% posw = perform_spherial_planar_sampling(pos_sphere, type)
%
%   'type' can be 'area' or 'gnomonic'.
%
%   This is used to produced spherical geometry images.
%   The sampling is first projected onto an octahedron and then unfolded on
%   a square.
%
%   Copyright (c) 2004 Gabriel Peyre

if nargin<2
    sampling_type = 'area';
end

% all 3-tuple of {-1,1}
a = [-1 1];
[X,Y,Z] = meshgrid(a,a,a);

anchor_2d = {};
anchor_3d = {};
for i=1:8
   	x = X(i); y = Y(i); z = Z(i);
   	a2d = []; a3d = [];
	if z>0
		a2d = [a2d, [0;0]];
    else
		a2d = [a2d, [x;y]];
    end
    a2d = [a2d, [x;0]];
    a2d = [a2d, [0;y]];
    anchor_2d{i} = a2d;    
    a3d = [a3d, [0;0;z]];
    a3d = [a3d, [x;0;0]];
    a3d = [a3d, [0;y;0]];
    anchor_3d{i} = a3d;
end

pos = pos_sphere;
n = size(pos, 2);
posw = zeros(2,n);
for s = 1:8
    x = X(s); y = Y(s); z = Z(s);
    anc2d = anchor_2d{s};
    anc3d = anchor_3d{s};
    I = find( signe(pos(1,:))==x & signe(pos(2,:))==y & signe(pos(3,:))==z );
    posI = pos(:,I);
    nI = length(I);
    if strcmp(sampling_type, 'area')
        % find the area of the 3 small triangles
        p1 = repmat(anc3d(:,1), 1, nI);
        p2 = repmat(anc3d(:,2), 1, nI);
        p3 = repmat(anc3d(:,3), 1, nI);
        a1 = compute_spherical_area( posI, p2, p3 );
        a2 = compute_spherical_area( posI, p1, p3 );
        a3 = compute_spherical_area( posI, p1, p2 );
        % barycentric coordinates
        a = a1+a2+a3;
        % aa = compute_spherical_area( p1, p2, p3 );
        a1 = a1./a; a2 = a2./a; a3 = a3./a;
    elseif strcmp(sampling_type, 'gnomonic')
        % we are searching for a point y=b*x (projection on the triangle)
        % such that :   a1*anc3d(:,1)+a2*anc3d(:,2)+a3*anc3d(:,3)-b*x=0
        %               a1+a2+a3=1
        a1 = zeros(1,nI); a2 = a1; a3 = a1;
        for i=1:nI
            x = posI(:,i);
            M = [1 1 1 0; anc3d(:,1), anc3d(:,2), anc3d(:,3), x];
            a = M\[1;0;0;0];
            a1(i) = a(1);
            a2(i) = a(2);
            a3(i) = a(3);
        end
    else
        error('Unknown projection method.');
    end
    posw(:,I) = anc2d(:,1)*a1 + anc2d(:,2)*a2 + anc2d(:,3)*a3;
end


function y = signe(x)

y = double(x>=0)*2-1;

function A = compute_spherical_area( p1, p2, p3 )

% length of the sides of the triangles :
% cos(a)=p1*p2
a  = acos( dotp(p2,p3) );
b  = acos( dotp(p1,p3) );
c  = acos( dotp(p1,p2) );
s  = (a+b+c)/2;
% use L'Huilier's Theorem
% tand(E/4)^2 = tan(s/2).*tan( (s-a)/2 ).*tan( (s-b)/2 ).*tan( (s-c)/2 )
E = tan(s/2).*tan( (s-a)/2 ).*tan( (s-b)/2 ).*tan( (s-c)/2 );
A = 4*atan( sqrt( E ) );
A = real(A);

function d = dotp(x,y)
d = x(1,:).*y(1,:) + x(2,:).*y(2,:) + x(3,:).*y(3,:);



