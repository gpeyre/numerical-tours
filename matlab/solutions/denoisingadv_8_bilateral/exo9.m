glist = [1 1.2 2 4];
clf;
for i=1:4
    gamma = glist(i);
    fr = f1 + gamma*r;
    imageplot(clamp(fr), ['\gamma=' num2str(gamma)], 2,2,i);
end
