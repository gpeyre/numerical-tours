[tmp,i] = min(E); 
lambda = lambda_list(i);
hf = h(f,lambda);
clf;
imageplot(clamp(hf), strcat(['SNR=' num2str(snr(f0,hf),3) 'dB']));
