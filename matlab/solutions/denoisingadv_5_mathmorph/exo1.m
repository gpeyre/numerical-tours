clf;i = 1;
for w=[1 1.5 3 7]
    imageplot(strel(w), ['w=' num2str(w)], 1,4, i);
    i = i+1;
end
