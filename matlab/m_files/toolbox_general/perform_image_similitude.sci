function M1 = perform_image_similitude(M,mapping_type, u,u1,v,v1, w,w1)

// perform_image_similitude
//
//   M1 = perform_image_similitude(M,mapping_type,u,u1,v,v1,w,w1);
//
// Compute the affine similitude that map u to u1
// and v to v1, and then resample the image M.
// p and p1 are assumed to be in [0,1]²
//
//   If mapping_type=='similitude', compute a true similitude
//       T(x,y) = [a -b] * [x] + [c]
//                [b  a]   [y]   [d]
//   Solve the equations T(u)=u1 and T(v)=v1.
//
//   If mapping_type=='similitude', compute a true similitude
//       T(x,y) = [a  b] * [x] + [e]
//                [c  d]   [y]   [f]
//   Solve the equations T(u)=u1 and T(v)=v1 and T(w)=w1.
//
//   Copyright (c) 2006 Gabriel Peyré



////////////////////////////////////////////////////////////////////////
////// compute T //////
//   T(x,y) = Q * [x;y] + t

if strcmp(mapping_type, 'similitude')
    
    // the matrix of the linear system
    A = [u(1) -u(2) 1 0; u(2)  u(1) 0 1; v(1) -v(2) 1 0; v(2)  v(1) 0 1];
    // the right hand size
    rhs = [u1(:); v1(:)];
    // solve
    z = A \ rhs;
    // the similitude
    Q = [z(1) -z(2); z(2) z(1)];
    // the translation
    t = [z(3); z(4)];
    
elseif strcmp(mapping_type, 'affine')
        
    // the matrix of the linear system
    A = [u(1) u(2) 0    0    1 0; 0    0    u(1) u(2) 0 1; v(1) v(2) 0    0    1 0; 0    0    v(1) v(2) 0 1; w(1) w(2) 0    0    1 0; 0    0    w(1) w(2) 0 1];
    // the right hand size
    rhs = [u1(:); v1(:); w1(:)];
    // solve
    z = A \ rhs;
    // the similitude
    Q = [z(1) z(2); z(3) z(4)];
    // the translation
    t = [z(5); z(6)];
 
else
    error('Unknown mapping');
end


////////////////////////////////////////////////////////////////////////
////// perform resampling //////

// original grid in the warped domain
n = size(M,1);
x = linspace(0,1,n);
[X,Y] = meshgrid(x,x);
// inverse warping P1=T^-1(P)
P = [X(:)'; Y(:)']; // position of the sampling
P(1,:) = P(1,:) - t(1); // substract translation
P(2,:) = P(2,:) - t(2);
P1 = (Q^(-1))*P; // undo similitude
// reshape the results
X1 = reshape(P1(1,:),[n,n]);
Y1 = reshape(P1(2,:),[n,n]);


M1 = interp2(X,Y,M',X1,Y1)';

endfunction