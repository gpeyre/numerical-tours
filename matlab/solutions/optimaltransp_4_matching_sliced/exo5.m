tlist = linspace(0,1,6);
clf;
for i=1:length(tlist)
    t = tlist(i);
    ft = (1-t)*f + t*f1;
    subplot(2,length(tlist)/2,i);
    imageplot( func(hist2d(ft(1:2,:),Q)) );
    title(['t=' num2str(t,2)]);
end
