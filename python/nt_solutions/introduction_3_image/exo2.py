k = round(.8*n)
k = round(k/2)*2 # even number
Mf = fft2(M)
r = int(k/2)
Mf[m-r+2:m+r, m-r+2:m+r] = 0
Mh = real( ifft2(Mf) )
# display
clf
imageplot( M[m-20:m+20,m-20:m+20], 'Image', [1,2,1])
imageplot( Mh[m-20:m+20,m-20:m+20], 'Low pass filtered', [1,2,2])
