m = round(.1*n/3)*3;
pvertex = vertex*U;
% non linear 
pvertexN = perform_thresholding(pvertex,m,'largest');
vertexN = pvertexN*U';
% linear
pvertexL = pvertex;
pvertexL(:,m/3+1:n) = 0;
vertexL = pvertexL*U';
% display
clf;
subplot(1,2,1);
plot_mesh(vertexL,faces);
subplot(1,2,2);
plot_mesh(vertexN,faces);
disp(['Linear:     SNR=' num2str(snr(vertex,vertexL),3) 'dB']);
disp(['Non-linear: SNR=' num2str(snr(vertex,vertexN),3) 'dB']);
