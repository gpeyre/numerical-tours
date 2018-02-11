T = 3*sigma
w = 4
Mspin = np.zeros([n,n,n])

for i in range(w**3):
    # shift the image
    MnoisyC = circshift(Mnoisy, [dX[i],dY[i],dZ[i]])
    # denoise 
    MW = perform_haar_transf(MnoisyC, 1, + 1)
    MWT = perform_thresholding(MW, T, "hard")
    M1 = perform_haar_transf(MWT, 1, -1)
    # back
    M1 = circshift(M1, [-dX[i],-dY[i],-dZ[i]])
    # average the result
    Mspin = Mspin*(i)/(i+1) + M1/(i+1)