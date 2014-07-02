function J = plot_curvelet(MW, options)

% plot_curvelet - display curvelets coefficients
%
%   J = plot_curvelet(MW);
%
%   Based on curvelab.

%generate curvelet image (a complex array)
I = fdct_wrapping_dispcoef(MW);

% remove bckgd
U = (I==.5);

J = ones(size(I)+2);
JU = ones(size(I)+2);

J(2:end-1,2:end-1) = I;
JU(2:end-1,2:end-1) = U;

J(J==.5) = 1;
J = sign(J) .* abs(J).^.6;


if nargout==0
    hold on;
    imageplot(J);
    [c,h] = contour( JU, [.5 .5], 'r');
    set(h, 'LineWidth', 2);
end


function img = fdct_wrapping_dispcoef(C)

% fdct_wrapping_dispcoef - returns an image containing all the curvelet coefficients
%
% Inputs
%     C         Curvelet coefficients
%
% Outputs
%     img       Image containing all the curvelet coefficients. The coefficents are rescaled so that
%       the largest coefficent in each subband has unit norm.
%

[m,n] = size(C{end}{1});
nbscales = floor(log2(min(m,n)))-3;

img = rescale(C{1}{1},-1,1);
% img = img/max(max(abs(img))); %normalize
for sc=2:nbscales-1
    nd = length(C{sc})/4;
    wcnt = 0;
    
    ONE = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
        ONE = [ONE, fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w})];
    end
    wcnt = wcnt+nd;
    
    TWO = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
        TWO = [TWO; fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w})];
    end
    wcnt = wcnt+nd;
    
    THREE = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
        THREE = [fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w}), THREE];
    end
    wcnt = wcnt+nd;
    
    FOUR = [];
    [u,v] = size(C{sc}{wcnt+1});
    for w=1:nd
        FOUR = [fdct_wrapping_dispcoef_expand(u,v,C{sc}{wcnt+w}); FOUR];
    end
    wcnt = wcnt+nd;
    
    [p,q] = size(img);
    [a,b] = size(ONE);
    [g,h] = size(TWO);
    m = 2*a+g;    n = 2*h+b; %size of new image
    scale = max(max( max(max(abs(ONE))),max(max(abs(TWO))) ), max(max(max(abs(THREE))), max(max(abs(FOUR))) )); %scaling factor
    
    new = 0.5 * ones(m,n);%background value
    new(a+1:a+g,1:h) = FOUR/scale;
    new(a+g+1:2*a+g,h+1:h+b) = THREE/scale;
    new(a+1:a+g,h+b+1:2*h+b) = TWO/scale;
    new(1:a,h+1:h+b) = ONE/scale;%normalize
    
    dx = floor((g-p)/2);    dy = floor((b-q)/2);
    
    new(a+1+dx:a+p+dx,h+1+dy:h+q+dy) = img;
    
    img = new;
end

function A = fdct_wrapping_dispcoef_expand(u,v,B)
A = zeros(u,v);
[p,q] = size(B);
A(1:p,1:q) = B;



