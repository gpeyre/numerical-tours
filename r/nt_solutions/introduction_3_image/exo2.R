k = round(.8 * n)
k = round(k/2) * 2 # even number
Mf = fft(M[,])
Mf[c((n/2) - (k/2) + 2) : c((n/2) + (k/2)), c((n/2) - (k/2) + 2) : c((n/2) + (k/2))] = 0
Mh = Re(fft(Mf, inverse=TRUE) / n)
# display
imageplot(M[c(((n/2) - 25):((n/2) + 25)), c(((n/2) - 25):((n/2) + 25))], 'Image', c(1,2,1))
imageplot( Mh[c(((n/2) - 25):((n/2) + 25)), c(((n/2) - 25):((n/2) + 25))], 'Low pass filtered', c(1,2,2))