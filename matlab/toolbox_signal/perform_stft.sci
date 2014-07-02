function y = perform_stft(x,w,q, options)

// perform_stft - compute a local Fourier transform
//
// Forward transform:
//   MF = perform_stft(M,w,q, options);
// Backward transform:
//   M  = perform_stft(MF,w,q, options);
//
//   w is the width of the window used to perform local computation.
//   q is the spacing betwen each window.
//
//   MF(:,i) contains the spectrum around point (i-1)*q
//
//   A typical use, for an redundancy of 2 could be w=2*q+1
//
//   options.bound can be either 'per' or 'sym'
//
//   options.normalization can be set to
//       'tightframe': tight frame transform, with energy conservation.
//       'unit': unit norm basis vectors, usefull to do thresholding
//
//   Copyright (c) 2006 Gabriel Peyre

options.null = 0;
if size(x,1)==1 | size(x,2)==1
    x = x(:); direction = 1;
    n = length(x);
else
    direction = -1;
    n = getoptions(options, 'n', 1, 1);
end


bound = getoptions(options, 'bound', 'per');
transform_type = getoptions(options, 'transform_type', 'fourier');
normalization = getoptions(options, 'normalization', 'tightframe');
window_type = getoptions(options, 'window_type', 'sin');
eta = getoptions(options, 'eta', 1);

// perform sampling
X = 1:q:n+1;
p = length(X);

if pmodulo(w,2)==1
// w = ceil((w-1)/2)*2+1;
    w1 = (w-1)/2;
    dX = (-w1:w1)';
else
    dX = (-w/2+1:w/2)';
end

X1 = repmat(X, [w 1]) + repmat(dX, [1 p]);
if strcmp(bound,'sym')
        X1(X1<1) = 1-X1(X1<1);
        X1(X1>n) = 2*n+1-X1(X1>n);
elseif strcmp(bound,'per')
        X1 = pmodulo(X1-1,n)+1;
else
    error('Unkown boundary conditions.');
end
I = X1;

// build a weight function
if strcmp(window_type,'sin') | strcmp(window_type,'hanning')
//        t = linspace(-%pi,%pi,w);
//        W = cos(t(:))+1;
        W = .5 *(1 - cos( 2*%pi*(0:w-1)'/(w-1) ));
elseif strcmp(window_type,'constant')
        W = ones(w,1);
else
        error('Unkwnown winow.');
end


//// renormalize the windows
weight = zeros(n,1);
for i=1:p
    weight(I(:,i)) = weight(I(:,i)) + W.^2;
end
weight = sqrt(weight);
Weight = repmat(W, [1 p]);
for i=1:p
    Weight(:,i) = Weight(:,i) ./ weight(I(:,i));
end


if strcmp(normalization, 'unit')
    if strcmp(transform_type, 'fourier')
        // for Fourier it is easy
        Renorm = sqrt( sum( Weight.^2, 1 ) )/w;
    else
        error('Not yet implemented');
        // for DCT it is less easy ...
        // take a typical window in the middle of the image
        weight = Weight(:,:,round(end/2),round(end/2));
        // compute diracs
        [X,Y,fX,fY] = ndgrid(0:w-1,0:w-1,0:w-1,0:w-1);
        A = 2 * cos( pi/w * ( X+1/2 ).*fX ) .* cos( pi/w * ( Y+1/2 ).*fY ) / w;
        A(:,:,1,:) = A(:,:,1,:) / sqrt(2); // scale zero frequency
        A(:,:,:,1) = A(:,:,:,1) / sqrt(2); 
        A = A .* repmat( weight, [1 1 w w] );
        Renorm = sqrt( sum( sum( abs(A).^2, 1 ),2  ) );
    end
end
    
//// compute the transform
if direction==1
    y = zeros(eta*w,p);
    if pmodulo(w,2)==1
        m = (eta*w+1)/2; w1 = (w-1)/2;
        sel = m-w1:m+w1;
    else
        m = (eta*w)/2+1; w1 = w/2;
        sel = m-w1:m+w1-1;
    end
    y(sel,:) = reshape(x(I),size(I)) .* Weight;
    // perform the transform
    y = my_transform( y, +1, transform_type );
    // renormalize if necessary
    if strcmp(normalization, 'unit')
        y = y ./ repmat( Renorm, [1 p] );
    end
else
    if strcmp(normalization, 'unit')
        x = x .* repmat( Renorm, [1 p] );
    end
    x = my_transform( x, -1, transform_type );
    x = real( x.*Weight );
    y = zeros(n,1);
    for i=1:p
        y(I(:,i)) = y(I(:,i)) + x(:,i);
    end
end


endfunction


////////////////////////////////////////////
////////////////////////////////////////////
////////////////////////////////////////////
function y = my_transform(x,direction,transform_type)

// my_transform - perform either FFT or DCT with energy conservation.
//   Works on array of size (w,w,a,b) on the 2 first dimensions.

w = size(x,1);
if strcmp(transform_type, 'fourier')
    // normalize for energy conservation
    if direction==1
        y = myfft(x) / sqrt(w);
    else
        y = myifft( x*sqrt(w) );
    end
elseif strcmp(transform_type, 'dct')
    for i=1:size(x,2)
        y(:,i) = perform_dct_transform(x(:,i),direction);
    end
else
    error('Unknown transform');
end

endfunction

////////////////////////////////////////////
////////////////////////////////////////////
function y = myfft(x)

y = x;
for i=1:size(x,2)
    y(:,i) = fft(x(:,i));
end

endfunction
////////////////////////////////////////////
////////////////////////////////////////////
function y = myifft(x)

y = x;
for i=1:size(x,2)
    y(:,i) = ifft(x(:,i));
end

endfunction
