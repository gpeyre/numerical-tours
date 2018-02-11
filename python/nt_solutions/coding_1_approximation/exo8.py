fL = np.zeros([n, n])

for i in range(1,n//w+1):
    for j in range(1,n//w+1):
        fL[(i-1)*w: i*w, (j-1)*w: j*w] = dct2(f[(i-1)*w: i*w, (j-1)*w: j*w])