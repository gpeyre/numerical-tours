clf
bar(arange(1,p+1), abs(wSparse))
bar(arange(1,p+1), -abs(wRidge))
legend(('Lasso', 'Ridge'))
