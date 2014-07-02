function plot_geometry_image(M, display_cut, display_octa)

% plot_geometry_image - plot a geometry image
%
%   plot_geometry_image(M, display_cut, display_octa);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    display_cut = 0;
end
if nargin<3
    display_octa = 0;
end

M = rescale(M);
if ndims(M)==2
    M = M(:, end:-1:1);
end

n = size(M,1);
hold on;
if ndims(M)==2
    surf(1:n, 1:n, M*20);
    axis([1 n 1 n 0 40]);
else
    C = ones(size(M(:,:,1)))*0.8;
    surf(M(:,:,1), M(:,:,2), M(:,:,3), C );
    axis tight;
end

% display cut
if display_cut
    str = 'r'; w = 2;
    h = plot3( M(1:end/2,1,1), M(1:end/2,1,2), M(1:end/2,1,3), str );
    set(h,'LineWidth',w);
    h = plot3( M(1:end/2,end,1), M(1:end/2,end,2), M(1:end/2,end,3), str );
    set(h,'LineWidth',w);
    h = plot3( M(1,1:end/2,1), M(1,1:end/2,2), M(1,1:end/2,3), str );
    set(h,'LineWidth',w);
    h = plot3( M(end,1:end/2,1), M(end,1:end/2,2), M(end,1:end/2,3), str );
    set(h,'LineWidth',w);
end
if display_octa
    str = 'r'; w = 2;
    n = size(M,1);
    t = 0:1/(n-1):1;
    x = -1:2/(n-1):1;
    [Y,X] = meshgrid(x,x);
    
    plot_cut(0*x,x, X,Y, M, str);
    plot_cut(x,0*x, X,Y, M, str);
    plot_cut(t,1-t, X,Y, M, str);
    plot_cut(t,t-1, X,Y, M, str);
    plot_cut(t-1,t, X,Y, M, str);
    plot_cut(t-1,-t, X,Y, M, str);
end
hold off;

shading interp;
lighting phong;
% lighting none;
camlight local; 
camlight infinit; 
material dull;
material shiny;
camproj('perspective');
axis square;
axis off;
% colormap gray(256);
cameramenu;
view(-170, 55);
axis tight;


%%%
function plot_cut(x,y, X,Y, M, str)
w = 2;
a = interp2(Y,X,M(:,:,1),x,y);
b = interp2(Y,X,M(:,:,2),x,y);
c = interp2(Y,X,M(:,:,3),x,y);
h = plot3(a,b,c,str);
set(h,'LineWidth',w);