function [H,x,xc] = compute_conditional_histogram(M,C, options)

% compute_conditional_histogram - compute conditional histograms
%
% [H,x,xc] = compute_conditional_histogram(M,C, options);
%
%   H is the conditional histograms:
%       H(i,j) = P( M=x(i,j) | C=xc(i) )
%
%   options.p is the number of bins used for M
%   options.pc is the number of bins used for C
%
%   If M (resp. C) are integer then p and pc are fit to the range of M (resp. C).
%   If M (resp. C) contains negative values, then x (resp. xc) are
%       symmetric around zero.
%   Otherwise x and xc spans evenly the range of M and C.
%
%   Copyright (c) 2005 Gabriel Peyr?


options.null = 0;
vmin  = min( M(:) );
vmax  = max( M(:) );
vminc = min( C(:) );
vmaxc = max( C(:) );

flagint  = sum( abs(M(:)-round(M(:))) )==0;
flagintc = sum( abs(C(:)-round(C(:))) )==0;
flagneg  = sum(M(:)<0)>0;
flagnegc = sum(M(:)<0)>0;

p = 21;
if not(isfield(options, 'p'))
    if flagint & flagneg
        vmax = max(abs(M(:)));
        vmin = -vmax;
    end
    if flagint
        p = vmax-vmin+1;
    end
end
pc = 21;
if not(isfield(options, 'pc'))
    if flagintc & flagnegc
        vmaxc = max(abs(C(:)));
        vminc = -vmaxc;
    end
    if flagintc
        pc = vmaxc-vminc+1;
    end
end
if isfield(options,'normalize')
    normalize = options.normalize;
else
    normalize = 1;
end

sc = (vmaxc-vminc)/(pc-1);
s = (vmax-vmin)/(p-1);
xc = vminc:sc:vmaxc;
x = vmin:s:vmax;

H = zeros(p,pc);

for ic = 1:length(xc)
    I = find( abs(C-xc(ic))<=sc/2 );
    for i = 1:length(x)
        H(i,ic) = sum( M(I)==x(i) );
    end 
end
if normalize
    A = repmat( sum(H),[p 1] );
    A(A<eps) = eps;
    H = H ./ A;
end