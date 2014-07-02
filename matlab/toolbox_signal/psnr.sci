function p = psnr(x,y, vmax)

// psnr - compute the Peack Signal to Noise Ratio
//
//   p = psnr(x,y,vmax);
//
//   defined by :
//       p = 10*log10( vmax^2 / |x-y|^2 )
//   |x-y|^2 = mean( (x(:)-y(:)).^2 )
//   if vmax is ommited, then
//       vmax = max(max(x(:)),max(y(:)))
//
//   Copyright (c) 2004 Gabriel Peyre

if argn(2)<3
    m1 = max( abs(x(:)) );
    m2 = max( abs(y(:)) );
    vmax = max(m1,m2);
end

d = mean( (x(:)-y(:)).^2 );

p = 10*log10( vmax^2/d );

endfunction