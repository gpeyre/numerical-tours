function err = compute_l2err(MF, donormalize, isnonlinear)

// compute_l2err - non linear approximation error
//
//   err = compute_l2err(MF,donormalize);
//
// err(m) is the square L2 non linear
// approximation error.
//
//   Copyright (c) 2008 Gabriel Peyre

if argn(2)<2
    donormalize = 0;
end

v = sort(abs(MF(:)).^2);
if v($)<v(1)
    v = reverse(v); // be sure it is in reverse order
end
err = reverse( cumsum(v) ); 
if donormalize
    err = err/err(1);
end

endfunction