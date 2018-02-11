w = random.standard_normal([n0,n0])/n0
w = w-np.mean(w) + 1/n0**2

w_transf = pyl.fft2(w)
w_transf_extended = np.zeros([n0,n0,np.shape(f0)[2]])
for i in range(np.shape(f0)[2]):
    w_transf_extended[:,:,i] = w_transf
    
f = np.real(pyl.ifft2(pyl.fft2(f1, axes=(0,1))*w_transf_extended, axes=(0,1)))

plt.figure(figsize = (10,10))

u = np.ones([n0,n0,3])
u[n0//2 - n//2 :n0//2 + n//2, n0//2 - n//2:n0//2 + n//2, :] = f0
u[0,0,:] = 0
u[1,0,:] = 1
imageplot(clamp(u), 'Input', [1,2,1])

u = f
u[0,0,:] = 0
u[1,0,:] = 1
imageplot(clamp(u), 'Synthesized', [1,2,2])