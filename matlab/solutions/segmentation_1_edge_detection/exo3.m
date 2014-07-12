slist = [1 2 4 6];
clf;
for i=1:length(slist)
    sigma = slist(i);
    d = sqrt( sum(nabla(  blur(f0,sigma)  ).^2,3) );
    t = max(d(:)) * 1/5;
    imageplot( double(d>t), ['\sigma=' num2str(sigma)], 2,2,i );
end
