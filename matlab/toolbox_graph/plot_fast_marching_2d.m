function plot_fast_marching_2d(W,S,path,start_points,end_points, options)

% plot_fast_marching_2d - plot the result of the fast marching.
%
%   plot_fast_marching_2d(W,S,path,start_points,end_points, options);
%
%   Copyright (c) 2004 Gabriel Peyré


n = size(W,1);

options.null = 0;

path_width = getoptions(options, 'path_width', 3);
point_size = getoptions(options, 'point_size', 12);


if isfield(options, 'start_point_style')
    start_point_style = options.start_point_style;
else
    start_point_style = 'rx';
end
if isfield(options, 'end_point_style')
    end_point_style = options.end_point_style;
else
    end_point_style = 'k*';
end

if ~iscell(end_point_style)
    end_point_style = {end_point_style, 'k.'};
end
if ~iscell(point_size)
    point_size = {point_size, 20};
end

if size(W,3)==1
    Z = repmat(W',[1,1,3]);
else
    Z = permute(W,[2 1 3]);
end
S1 = repmat(S',[1,1,3]);
S1(:,:,1) = Inf;
open_list = find( S1==0 );
close_list = find( S1==-1 );

set(gca,'box','on');

hold on;
% Z(open_list) = 0;
Z(close_list) = 0;
imagesc( rescale(Z) );
axis image;
axis off;

c_list = perform_curve_extraction(S,0);
for k=1:length(c_list);
    c_list{k} = c_list{k}*(n-1)+1;
end
for i=1:length(c_list)
    plot(c_list{i}(1,:), c_list{i}(2,:), 'r-', 'LineWidth', 2);
end

if nargin>=3
    if ~iscell(path)
        path = {path};
    end
    for i=1:length(path)
        p = path{i};
        if size(p,1)>size(p,2)
           p = p'; 
        end
        plot(p(1,:), p(2,:), '-', 'LineWidth', path_width);
    end
end

% axis square;
if nargin>=4 & ~isempty(start_points)
    if size(start_points,1)~=2
        start_points = start_points';
    end
    if size(start_points,1)~=2
        error('start_points must be of size 2xk');
    end
    for i=1:size(start_points,2)
	    plot(start_points(1,i), start_points(2,i), start_point_style, 'MarkerSize', point_size{1});
    end
end
if nargin>=5 & ~isempty(end_points)
    if size(end_points,1)~=2
        end_points = end_points';
    end
    if size(end_points,1)~=2
        error('end_points must be of size 2xk');
    end
    for i=1:size(end_points,2)
        plot(end_points(1,i), end_points(2,i), end_point_style{1 + (i>1)}, 'MarkerSize', point_size{1 + (i>1)});
    end
end
hold off;