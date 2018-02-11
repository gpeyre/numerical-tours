k = round(.8*n)
k = round(k/2)*2 # even number
Mf = plan_fft(M)*M
Mf[Int(n/2 - k/2) + 3 : Int(n/2 + k/2), Int(n/2 - k/2) + 3 : Int(n/2 + k/2)] = 0
Mh = real( plan_ifft(Mf)*Mf )
# display
clf
imageplot( M[Int(n/2 - 19) : Int(n/2 + 20), Int(n/2 - 19) : Int(n/2 + 20)], "Image", [1, 2, 1])
imageplot( Mh[Int(n/2 - 19) : Int(n/2 + 20), Int(n/2 - 19) : Int(n/2 + 20)], "Low pass filtered", [1, 2, 2])
