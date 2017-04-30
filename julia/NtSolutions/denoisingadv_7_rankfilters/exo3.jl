figure(figsize = (5,5))
openingclosing = f -> closing(opening(f))
imageplot(openingclosing(f))
