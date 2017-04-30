figure(figsize = (14, 7))

f1 = f
for i in 1 : 4
    f1 = opening(f1)
    imageplot(f1, @sprintf("Iteration %i" ,i), [2, 4, i])
end

f1 = f
for i in 1 : 4
    f1 = closing(f1)
    imageplot(f1, @sprintf("Iteration %i", i), [2, 4, i + 4])
end
