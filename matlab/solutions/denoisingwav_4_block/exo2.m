tlist = linspace(.3,.9,4);
clf;
for i=1:length(tlist)
    T = tlist(i);
    imageplot( clamp(ThreshBlock(f,T)), ['T=' num2str(T,2)], 2,2,i );
end
