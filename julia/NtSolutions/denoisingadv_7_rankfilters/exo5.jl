## Insert your code here.

figure(figsize = (10, 7))

f1 = f
for i in 1 : 6
    f1 = medfilt(f1)
    imageplot(f1, @sprintf("Iteration %i", i), [2, 3, i])
end
