from nt_toolbox.perform_haar_transf import *
from nt_toolbox.perform_thresholding import *

MW = perform_haar_transf(Mnoisy, 1, +1)
Tlist = np.linspace(1, 4, 20)*sigma
err_hard = []
err_soft = []

for i in range(len(Tlist)):
    MWT = perform_thresholding(MW, Tlist[i], 'hard')
    M1 = perform_haar_transf(MWT, 1, -1)
    err_hard = err_hard + [snr(M, M1)]
    MWT = perform_thresholding(MW, Tlist[i], 'soft')
    M1 = perform_haar_transf(MWT, 1, -1)
    err_soft = err_soft + [snr(M, M1)]
    if i > 1 and err_soft[i] > np.max(err_soft[:i]):
        Mwav = M1

plt.figure(figsize=(7,5))

plt.plot(Tlist/ sigma, err_hard, '.-', label="hard", color="b")
plt.plot(Tlist/ sigma, err_soft, '.-', label="soft", color="r")
plt.xlabel("$T/\sigma$")
plt.ylabel("SNR")
plt.ylim(np.min(err_hard),np.max(err_soft))
plt.legend(loc='upper right')

plt.show()