f1 = f;
clf;
for i=1:6
    f1 = medfilt(f1);
    imageplot(f1, ['iteration ' num2str(i)], 2,3, i);
end
