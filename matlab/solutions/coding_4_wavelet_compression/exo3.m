Tlist = linspace(0.01,0.3,60);
for i=1:length(Tlist)
    err(i) = norm( M-perform_thresholding(M, Tlist(i), 'quantize'), 'fro');
end 
errW = norm(M-MQ,'fro');
[tmp,i] = min( abs(err-errW) );
T0 = Tlist(i);
MQ0 = perform_thresholding(M, T0, 'quantize');
disp(['Spatial quantization step T0=' num2str(T0,2) '.']);
clf
imageplot( clamp(MQ), 'Wavelet quantized',1,2,1);
imageplot( clamp(MQ0), 'Spacial quantized',1,2,2);
