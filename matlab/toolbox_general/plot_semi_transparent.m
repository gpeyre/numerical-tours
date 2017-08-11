function plot_semi_transparent(x,y,col,alpha)

% plot_semi_transparent - plot a superposition of transparent plots
%
%    plot_semi_transparent(x,y,col,alpha);
%
% each y(:,i) is a signal to plot, must have same length as x.
%
%   Copyright (c) Gabriel Peyre

hold on;
for i=1:size(y,2);
    a = [x fliplr(x)];
    b = [y(:,i)', fliplr(y(:,i)')];
	patch(a,b,'r','EdgeColor',col, 'EdgeAlpha',alpha,'FaceColor','none','Linesmoothing','on');
end
axis tight;


end