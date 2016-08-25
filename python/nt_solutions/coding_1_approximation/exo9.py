from numpy import linalg

f1 = np.copy(fL)

for i in range(1,n//w+1):
    for j in range(1,n//w+1):
        f1[(i-1)*w: i*w, (j-1)*w: j*w] = idct2(f1[(i-1)*w: i*w, (j-1)*w: j*w])

print("Error |f-f1|/|f| =", linalg.norm(f-f1)/linalg.norm(f))