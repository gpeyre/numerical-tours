




mydisp <- function(x){ log10(abs(fftshift_4d(x)) + 1e-5) }


# Original
FFT_1 <- apply(f, c(3,4), fft)
FFT_1_arr <- array(0, dim(f))
for (i in 1:n){ for (j in 1:n){ FFT_1_arr[i,j,,] <- FFT_1[(i-1)*n + j,,] } }

DISP_1 <- mydisp(FFT_1_arr)

imageplot(apply(DISP_1, c(1,2), mean), 'Original', c(1,2,1))


# Periodic layer
FFT_2 <- apply(p, c(3,4), fft)
FFT_2_arr <- array(0, dim(p))
for (i in 1:n){ for (j in 1:n){ FFT_2_arr[i,j,,] <- FFT_2[(i-1)*n + j,,] } }

DISP_2 <- mydisp(FFT_2_arr)

imageplot(apply(DISP_2, c(1,2), mean), 'Periodic layer', c(1,2,2))