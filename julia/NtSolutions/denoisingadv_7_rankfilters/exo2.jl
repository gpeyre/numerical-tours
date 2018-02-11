figure(figsize = (5,5))
closingopening = f -> opening(closing(f))
imageplot(closingopening(f))
