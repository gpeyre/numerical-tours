function x = perform_wavelet_transf(x, Jmin, direction, options)

// perform_wavelet_transf - peform fast wavelet transform
//
//   y = perform_wavelet_transf(x, Jmin, direction, options);
//
//   Implement 1D and 2D symmetric wavelets with symmetric boundary treatements, using
//   a lifting implementation.
//
//   h = options.filter gives the coefficients of the lifting filter.
//   You can use h='linear' or h='7-9' to select automatically biorthogonal
//   transform with 2 and 4 vanishing moments.
//
//   You can set options.ti=1 to compute a translation invariant wavelet
//   transform.
//
//  For TI transform, scilab might run out of memory. To increase it, use the function extend_stack_size.
//
//   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
h = getoptions(options, 'filter', '7-9');
h = lower(h);
if isstr(h)
	// P/U/P/U/etc the last coefficient is scaling
   if strcmp(h,'linear') | strcmp(h, '5-3')
           h = [1/2 1/4 sqrt(2)];
   elseif strcmp(h, '9-7') | strcmp(h, '7-9')
           h = [1.586134342 -.05298011854 -.8829110762 .4435068522 1.149604398];
   else
           error('Unknown filter');
   end
end

// detect dimensionality 
// d = getoptions(options,'ndims',-1);
d = -1;
if isfield(options, 'ndims')
    ndims = options.ndims;
end
if d<0
    if size(x,1)==1 | size(x,2)==1
        d = 1;
    else
        d = 2;
    end
end

if d==1 & size(x,1)==1
    warning('For 1D transform, the vector should be a column vector.');
    x = x(:);
end

ti = getoptions(options, 'ti', 0);

separable = getoptions(options, 'separable', 0);
if d==2 & separable==1
    options.ti = 0;
    if ti==1
        warning('Separable does not works for translation invariant transform');
    end
    // perform a separable wavelet transform
    n = size(x,1);
    if direction==1
        for i=1:n
            x(:,i) = perform_wavelet_transf(x(:,i),Jmin,direction,options);
        end
        for i=1:n
            x(i,:) = perform_wavelet_transf(x(i,:)',Jmin,direction,options)';
        end
    else
        for i=1:n
            x(i,:) = perform_wavelet_transf(x(i,:)',Jmin,direction,options)';
        end
        for i=1:n
            x(:,i) = perform_wavelet_transf(x(:,i),Jmin,direction,options);
        end
    end
    return;
end



// number of lifting steps
n = size(x,1);
m = (length(h)-1)/2;
Jmax = log2(n)-1;
jlist = Jmax:-1:Jmin;
if direction==-1
    jlist = jlist($:-1:1);
end

if ti==0
    // subsampled
    for j=jlist
        if d==1
            x(1:2^(j+1),:) = lifting_step( x(1:2^(j+1),:), h, direction);
        else
            x(1:2^(j+1),1:2^(j+1)) = lifting_step( x(1:2^(j+1),1:2^(j+1)), h, direction);
            x(1:2^(j+1),1:2^(j+1)) = lifting_step( x(1:2^(j+1),1:2^(j+1))', h, direction)';
        end
    end
else
    // TI
    nJ = Jmax-Jmin+1;
    if direction==1 & d==1
        x = repmat(x,[1 1 nJ+1]);
    elseif direction==1 & d==2
        x = repmat(x,[1 1 3*nJ+1]);
    end
    for j=jlist
        dist = 2^(Jmax-j);
        if d==1
            if direction==1
                x(:,:,[1 j-Jmin+2]) = lifting_step_ti( x(:,:,1), h, direction, dist);
            else
                x(:,:,1) = lifting_step_ti( x(:,:,[1 j-Jmin+2]), h, direction, dist);                
            end
        else
            dj = 3*(j-Jmin);
            if direction==1
                x(:,:,[1 dj+2]) = lifting_step_ti( x(:,:,1), h, direction, dist);
                
                x(:,:,[1 dj+3]) = lifting_step_ti( x(:,:,1)', h, direction, dist);
                x(:,:,1) = x(:,:,1)'; x(:,:,dj+3) = x(:,:,dj+3)';
                
                x(:,:,[dj+[2 4]]) = lifting_step_ti( x(:,:,dj+2)', h, direction, dist);
                x(:,:,dj+2) = x(:,:,dj+2)';
                x(:,:,dj+4) = x(:,:,dj+4)';                
            else                
                x(:,:,dj+2) = x(:,:,dj+2)';
                x(:,:,dj+4) = x(:,:,dj+4)';
                x(:,:,dj+2) = lifting_step_ti( x(:,:,[dj+[2 4]]), h, direction, dist)';

                x(:,:,1) = x(:,:,1)'; 
                x(:,:,dj+3) = x(:,:,dj+3)';
                x(:,:,1) = lifting_step_ti( x(:,:,[1 dj+3]), h, direction, dist)';
                
                x(:,:,1) = lifting_step_ti( x(:,:,[1 dj+2]), h, direction, dist);
           end            
        end
    end
    if direction==-1
        x = x(:,:,1);
    end
end


endfunction

// --- //
function x = lifting_step(x, h, direction)

// number of lifting steps
m = (length(h)-1)/2;

if direction==1
    // split
    d = x(2:2:$,:);
    x = x(1:2:$,:);
    k = size(x,1);
    for i=1:m
        d = d - h(2*i-1) * ( x + x([2:k,k],:) );
        x = x + h(2*i ) * ( d + d([1,1:k-1],:) );
    end
    x = [x*h($);d/h($)];
else
    // retrieve detail coefs
    d = x($/2+1:$,:)*h($);
    x = x(1:$/2,:)/h($);
    k = size(x,1);
    for i=m:-1:1
        x = x - h(2*i  ) * ( d + d([1 1:k-1],:) );
        d = d + h(2*i-1) * ( x + x([2:k k],:) );
    end
    // merge
    x1 = zeros(size(x,1)*2,size(x,2));
    x1(1:2:$,:) = x;
    x1(2:2:$,:) = d;
    x = x1;
end

endfunction

// --- //
function x = lifting_step_ti(x, h, direction, dist)

// number of lifting steps
m = (length(h)-1)/2;
n = size(x,1);

s1 = (1:n)+dist;
s2 = (1:n)-dist;
// boundary conditions
s1(s1>n) = 2*n-s1(s1>n); s1(s1<1) = 2-s1(s1<1);
s2(s2>n) = 2*n-s2(s2>n); s2(s2<1) = 2-s2(s2<1);

// s1 = [2 1:n-1]; s2 = [2:n n-1];

if direction==1
    // split
    d = x;
    for i=1:m
        d = d - h(2*i-1) * ( x(s1,:) + x(s2,:) );
        x = x + h(2*i  ) * ( d(s1,:) + d(s2,:) );
    end
    x = cat_3d(x*h($),d/h($));
else
    // retrieve detail coefs
    d = x(:,:,2)*h($);
    x = x(:,:,1)/h($);
    for i=m:-1:1
        x = x - h(2*i  ) * ( d(s1,:) + d(s2,:) );
        d = d + h(2*i-1) * ( x(s1,:) + x(s2,:) );
    end
    // merge
    x = (x+d)/2;
end

endfunction


// --- //
function c = cat_3d(a,b)
// I do not know why, but cat is not working properly on scilab ...
c = zeros(size(a,1), size(a,2), 2);
c(:,:,1) = a;
c(:,:,2) = b;

endfunction
