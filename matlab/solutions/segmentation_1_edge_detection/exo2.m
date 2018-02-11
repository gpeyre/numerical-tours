t_list = max(d(:)) * [1/4 1/5 1/10 1/20];
clf;
for i=1:length(t_list)
    t = t_list(i);
    imageplot( double(d>t), ['t=' num2str(t, 2)] , 2,2,i);
end
