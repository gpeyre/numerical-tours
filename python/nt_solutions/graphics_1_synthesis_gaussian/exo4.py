plt.figure(figsize = (10,10))

for i in range(4):
    w = random.standard_normal([n,n])/n
    w = w-np.mean(w) + 1/n**2
    
    w_transf = pyl.fft2(w)
    w_transf_extended = np.zeros([n,n,np.shape(f0)[2]])
    for k in range(np.shape(f0)[2]):
        w_transf_extended[:,:,k] = w_transf
        
    f = np.real(pyl.ifft2(pyl.fft2(f0, axes =(0,1))*w_transf_extended, axes =(0,1)))

    u = f
    u[0,0] = 0
    u[1,0] = 1
    imageplot(clamp(u), 'Synthesized', [2,2,i+1])