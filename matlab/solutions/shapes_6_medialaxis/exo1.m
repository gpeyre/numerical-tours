% gradient
G = grad(Q);
G(G<-nbound/2) = G(G<-nbound/2) + nbound;
G(G>nbound/2) = G(G>nbound/2) - nbound;
% Compute the norm of the gadient.
G = sqrt(sum(G.^2,3));
% Remove the boundary to the skeletton.
M1 = perform_convolution(M,ones(3)/9)>.99;
G = G.*M1;
clf;
B = G>15;
A = M;
A(B==1) = 0;
imageplot(-A);
