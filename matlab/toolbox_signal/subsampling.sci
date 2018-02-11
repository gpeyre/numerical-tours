function y = subsampling(x,d,p)

// downsampling - subsampling along dimension d
//
//   y = upsampling(x,d,p);
//
//   default is p==2, d==1
//
//   Copyright (c) 2009 Gabriel Peyre

if argn(2)<3
    p = 2;
end
if argn(2)<2
    d = 1;
end

if d==1
        y = x(1:p:$,:,:);
elseif d==2
        y = x(:,1:p:$,:);
elseif d==3
        y = x(:,:,1:p:$);
else
        error('Not implemented');
end

endfunction