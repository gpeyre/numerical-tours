function B = display_shape_function(A,options)

% display_shape_function - display a function inside a shape.
%
%   B = display_shape_function(A,options);
%
%   options.cm gives colormap.
%   options.display_levelsets = 1 to show level sets.
%   options.pstart to give location of start point.
%
%   Copyright (c) Gabriel Peyre 2010

options.null = 0;

cm = getoptions(options, 'cm', jet(256));
display_levelsets = getoptions(options, 'display_levelsets', 0);
nbr_levelsets = getoptions(options, 'nbr_levelsets', 20);
pstart = getoptions(options, 'pstart', []);

A(isinf(A)) = 0;
I = find(A>0);
A(I) = rescale(A(I));

n = size(A,1);
B = ones(n,n,3);

m = size(cm,1);


for i=1:3
    c = cm(:,i);
    Bi = B(:,:,i);   
    Bi(I) = c( round(A(I)*(m-1))+1 );
    B(:,:,i) = Bi;
end

if nargout==0
	hold on;
    imageplot(B);
    if display_levelsets
        contour(rescale(A), nbr_levelsets, 'k');
        % set(gca, 'LineWidth', 2);
    end
    if not(isempty(pstart))
        h = plot(pstart(2), pstart(1), 'r.');
        set(h, 'MarkerSize', 30);
    end
    hold off;
    axis ij;
end