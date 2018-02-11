clf;
f1 = f;
for i=1:4
    f1 = opening(f1);
    imageplot(f1, ['iteration ' num2str(i)], 2,4, i);
end
f1 = f;
for i=1:4
    f1 = closing(f1);
    imageplot(f1, ['iteration ' num2str(i)], 2,4, 4+i);
end
