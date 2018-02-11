function imageplot(M,str, a,b,c)

% imageplot - diplay an image and a title
%
% Example of usages:
%   imageplot(M);
%   imageplot(M,title);
%   imageplot(M,title,1,2,1);   % to make subplot(1,2,1);
%
%   If you want to display several images:
%       imageplot({M1 M2}, {'title1', 'title2'});
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<2
    str = [];
end

nbdims = 2;
if size(M,1)==1 || size(M,2)==1
    nbdims = 1;
end

if iscell(M)
    q = length(M);
    if nargin<5
        c = 1;
        a = ceil(q/4);
        b = ceil(q/a);
    end
    if (c-1+q)>(a*b)
        warning('a and c parameters not large enough');
        a = ceil((c-1+q)/4);
        b = ceil((c-1+q)/a);
    end
    for i=1:q
        if iscell(str)
            str1 = str{i};
        else
            str1 = str;
        end
        imageplot(M{i},str1, a,b,c-1+i);
    end
    global axlist;
    if not(isempty(axlist))
        linkaxes(axlist, 'xy');
    end
    return;
end

if nargin==5
    global axlist;
    global imageplot_size;
    if c==1 || isempty(imageplot_size) || imageplot_size~=size(M,1)
        clear axlist; 
        global axlist; 
        axlist = [];
        imageplot_size = size(M,1);
    end
    axlist(end+1) = subplot(a,b,c);
end



if nbdims==1
    plot(M); axis tight;
else
    if size(M,3)==2
        M = cat(3,M, zeros(size(M,1),size(M,2)));
    end
    if size(M,3)==1
        colormap gray(256);
    else
        colormap jet(256);
    end
    imagesc(rescale(M)); axis image; axis off;
end
if not(isempty(str))
    title(str);
end


if nargin==5 && c==a*b
    linkaxes(axlist, 'xy');
end