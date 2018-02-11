
# Defining a stemplot function

stemplot = function(x, y, pch=16, linecol=1, clinecol=1,...)
{
if (missing(y))
{
    y = x
    x = 1:length(x)
 }
plot(x,y,pch=pch,...)
for (i in 1:length(x))
{
  lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
}
    lines(c(x[1] - 2, x[length(x)] + 2), c(0,0),col=clinecol)
}
