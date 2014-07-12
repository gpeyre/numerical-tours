clf;
i = 0;
for w = [1 1.5 2 4]
    i = i+1;
    imageplot(opening(M,w), ['w=' num2str(w)], 2,2,i);
end
