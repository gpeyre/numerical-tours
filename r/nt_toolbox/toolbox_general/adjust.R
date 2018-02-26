rotate = function(x) t(apply(x, 1, rev))
adjust = function(M){t(rotate(t(M)))}