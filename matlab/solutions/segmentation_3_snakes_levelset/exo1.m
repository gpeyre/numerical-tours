do_fast = 1;
% radius
r = n/3;
% center
c = n  - 10 - [r r];
% shape
phi2 = max( abs(X-c(1)), abs(Y-c(2)) ) - r;
