function A = circshift(A,s)

// circshift - circular shift an array
//
//  A = circshift(A,s);  
//
//  Copyright (c) 2009 Gabriel Peyre

if length(s)==2
    A = circshift(A,s(1));
    A = circshift(A',s(2))';
    return;
end

s = -s;

sel = 1:size(A,1);
sel = mod(sel+s-1,size(A,1))+1;
d = nb_dims(A);
if d==1
    A = A(sel);
elseif d==1
    A = A(sel,:);
elseif d==2
    A = A(sel,:,:);
elseif d==3
    A = A(sel,:,:);
elseif d==4
    A = A(sel,:,:,:);
end

endfunction