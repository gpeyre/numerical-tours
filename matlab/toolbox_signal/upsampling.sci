function y = upsampling(x,d,p)

// upsampling - add p zeros between samples along dimension d
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
    if ndims(x)==1
        y = zeros(p*size(x,1),1); 
    elseif ndims(x)==2
        y = zeros(p*size(x,1),size(x,2)); 
    else
        y = zeros(p*size(x,1),size(x,2),size(x,3)); 
    end
    y(1:p:size(y,1),:,:) = x;
elseif d==2
    if ndims(x)==3
        y = zeros(size(x,1),p*size(x,2),size(x,3)); 
    else
        y = zeros(size(x,1),p*size(x,2)); 
    end
    y(:,1:p:size(y,2),:) = x;
elseif d==3
    y = zeros(size(x,1),size(x,2),p*size(x,3)); 
    y(:,:,1:p:size(y,3)) = x;
else
    error('Not implemented');
end

endfunction