k = 100
t = np.linspace(0,.5,k)

plt.figure(figsize = (10,10))

plt.subplot(2, 1, 1)
plt.bar(t[:-1], np.histogram(f0,t)[0]*k/n**2, width = 1/200, color = "darkblue")
plt.title('Input')

plt.subplot(2, 1, 2)
plt.bar(t[:-1], np.histogram(f,t)[0]*k/n**2, width = 1/200, color = "darkblue")
plt.title('Synthesized')

plt.show()