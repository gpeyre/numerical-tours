%% Image Compression with Wavelets
% This numerical tour uses wavelets to perform image compression.
% We consider a simple model for compression, where we only
% estimate the number of bits of the compressed data, without
% really performing the actual entropic coding.

perform_toolbox_installation('signal', 'general');

%% Wavelet Domain Quantization
% Image compression is perfomed by first quantizing the wavelet coefficients of
% an image.

%%
% A scalar quantizer of step size |T| uses the function |floor|. It has a
% twice larger zero bins.

%%
% Create values evenly spaced for quantization.

v = linspace(-1,1, 2048);


%%
% Bin size for the quantization. The larger, the more agressive the
% compression.

T = .1;

%%
% For compression, we compute quantized integer values.

vI = floor(abs(v/T)).*sign(v);

%%
% For decompression, we compute de-quantized values from |vI|,
% which are chosen as the mid-point of each quantization bin.

vQ = sign(vI) .* (abs(vI)+.5) * T;

%%
% Display the quantization curve. 

clf;
subplot(1,2,1);
plot(v, vI);
axis('tight');
title(strcat(['Quantized integer values, T=' num2str(T)]));
subplot(1,2,2);
hold('on');
plot(v, vQ); 
plot(v, v, 'r--');
axis('equal'); axis('tight');
title('De-quantized real values');

%% Quantization and Approximation of Wavelet Coefficients
% Quantization of wavelet coefficients set to 0 those coefficients
% which are smaller than |T|, but it also modify the values of larger
% coeffiients. It thus creates an error that is slightly larger than 
% simply performing an approximation with hard thresholding 
% at |T|.

%%
% First we load an image.

n = 256;
M = rescale( load_image('lena', n) );

%%
% Compute its wavelet transform.

Jmin = 4;
MW = perform_wavelet_transf(M,Jmin, +1);

%EXO
%% Compute the coefficients |MWT| obtained by thresholding at
%% |T| the coefficients |MW|. Compute the coefficients |MWQ| obtained
%% by quantizing with bin size |T| the same coefficients.
%% Display them using the function |plot_wavelet|.
% thresholding approximation
MWT = perform_thresholding(MW, T, 'hard');
MWQ = perform_thresholding(MW, T, 'quantize');
% display
clf;
subplot(1,2,1);
plot_wavelet(MWT,Jmin);
title('Thresholded coefficients');
subplot(1,2,2);
plot_wavelet(MWT-MWQ,Jmin);
title('Thresholded - Quantized');
%EXO

%EXO
%% Compare the effect of quantizing at |T=.2| and thresholding at |T=.2|
%% the wavelet coefficients of an image.
% inverse transform
MT = perform_wavelet_transf(MWT,Jmin, -1);
MQ = perform_wavelet_transf(MWQ,Jmin, -1);
% error
eT = snr(M,MT);
eQ = snr(M,MQ);
% display
clf;
imageplot(MT, strcat(['Thresholding, SNR=' num2str(eT,2)]), 1,2,1)
imageplot(MT-MQ, strcat(['Thresholding - Approximating, SNR=+' num2str(eT-eQ,2)]), 1,2,2);
%EXO

%EXO
%% Compute a bin size |T0| to quantize the original |M| itself to obtained
%% |MQ0| such that |norm(M-MQ,'fro')| is as close as possible to the error
%% obtained with wavelet domain quantization.
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
%EXO

%% Entropy Coding the Wavelet Coefficients
% Actually store the quantized coefficients in a file, one need to compute
% a binary code from |MWI|. The length of this code is the number of bits 
% used by the compressor, which typically increases when |T| decays toward
% 0.

%%
% To reduce the number of bits, an entropic coder makes use of the
% statistical distribution of the quantized values.

%%
% First we quantize the coefficients.

MWI = floor(abs(MW/T)).*sign(MW);
MWQ = sign(MWI) .* (abs(MWI)+.5) * T;

%%
% Assuming that all the coefficients of |MWI| are drawn independently from the
% same distribution with histogram |h|, the minium bit per pixel achievable
% is the Entropy lower bound.

%%
% |E = -\sum_i \log2(h(i))*h(i)|

%%
% Huffman trees and more precisely block-Huffman tree codes get
% increasingly closer to this bound when the data size increases.
% Arithmetic coders also achieves very good results and are fast to
% compute.

%%
% Compute the nomalized histogram of the 
% quantized wavelet coefficients.

a = max(abs(MWI(:))); 
t = -a:a;
h = hist(MWI(:), t); h = h/sum(h);

%%
% Compute the histogram of the quantized pixels or the original image.

t0 = 0:1/T0;
MI = floor(abs(M/T0)); % quantized pixel values
h0 = hist(MI(:), t0); h0 = h0/sum(h0);

%%
% Display the histograms.

clf;
subplot(2,1,1);
bar(t0,h0); axis('tight');
title('Pixels');
subplot(2,1,2);
bar(t,h); axis([-5 5 0 max(h)])
title('Wavelets (zoom)');

%EXO
%% Compute the entropy lower bound for the quantized 
%% wavelet coefficients and for the quantized pixel values. 
%% Take care of |log2(0)| when |h(i)=0|.
h = h+1e-10; h = h/sum(h);
E = -sum( h.*log2(h) );
h0 = h0+1e-10; h0 = h0/sum(h0);
E0 = -sum( h0.*log2(h0) );
disp(['Pixels entropy:  ' num2str(E0,2)]);
disp(['Wavelet entropy: ' num2str(E,2)]);
%EXO

%EXO
%% Compute, for various threshold |T|, the number of bits per pixels |E(T)|
%% of the quantized wavelet coefficients,
%% and the wavelet decompression error |err(T)|, compute using SNR. 
%% Display the rate
%% distortion curve |err| as a function of |E|.
Tlist = linspace(.03,.6,20);
err = []; nbits = [];
for i=1:length(Tlist)
    T = Tlist(i);
    % quantize
    MWI = floor(abs(MW/T)).*sign(MW);
    MWQ = sign(MWI) .* (abs(MWI)+.5) * T;
    % inverse
    MQ = perform_wavelet_transf(MWQ,Jmin, -1);
    % error
    err(i) = snr(M,MQ);
    % bits
    nbits(i) = compute_entropy(MWI(:));
end
clf;
hh = plot(nbits,err); axis('tight');
set_label('bpp','SNR');
if using_matlab()
    set(hh,'LineWidth',2);
end
%EXO

%% Scale-by-scale Entropy Coding
% Wavelet coefficients of an image does not have the same distribution
% accross the scales. Taking this into account can further reduce the
% number of bits for coding.

%%
% Quantize the coeffients.

T = .1;
MWI = floor(abs(MW/T)).*sign(MW);

%%
% Extact the fine scale wavelet coefficients.

MWH = MWI(1:n/2,n/2+1:n);
MWV = MWI(n/2+1:n,1:n/2);
MWD = MWI(n/2+1:n,n/2+1:n);

%%
% Display.

clf;
imageplot(MWH,'Horizontal',1,3,1);
imageplot(MWV,'Vertical',1,3,2);
imageplot(MWD,'Diagonal',1,3,3);

%EXO
%% Extract the three fine scale wavelet coefficients (horizontal, vertical,
%% diagonal directions) and quantize them, for instance with |T=.1|.
%% Compute the entropy of the three sets together, and compute the entropy
%% of each set. 
Etot = compute_entropy([MWH MWV MWD]);
Ehor = compute_entropy(MWH);
Ever = compute_entropy(MWV);
Edia = compute_entropy(MWD);
disp([ 'Entropy, all:  ' num2str(Etot,3) ]);
disp([ 'Entropy, hor:  ' num2str(Ehor,3) ]);
disp([ 'Entropy, vert: ' num2str(Ever,3) ]);
disp([ 'Entropy, diag: ' num2str(Edia,3) ]);
%EXO

%EXO
%% Compare the number of bits needed to code all the wavelet coefficients
%% together, and the number of bits needed to code independantly each scale
%% of wavele coefficients for |Jmin=4<=j<=log2(n)-1| (and group together the
%% remaining coefficients for |j<Jmin|).
Esep = 0;
Jmax = log2(n)-1; Jmin = 4; 
for j = Jmax:-1:Jmin
    for q=1:3
        [selx,sely] = compute_quadsel(j,q);
        MWj = MWI(selx,sely);
        Esep = Esep + prod(size(MWj))*compute_entropy(MWj);
    end
end
Esep = Esep + prod(size(MWj))*compute_entropy(MWI(1:2^j,1:2^j));
Ewhole = compute_entropy(MWI);
disp(['nb.bis, whole:    ' num2str(Ewhole,3) ' bpp']);
disp(['nb.bis, separate: ' num2str(Esep/n^2,3) ' bpp']);
%EXO


%EXO
%% Compute the rate distortion curve obtained by coding the coefficient
%% separately through the scale, and compare with the rate distortion curve
%% obtained by coding the coefficients as a whole.
nbits1 = [];
for i=1:length(Tlist)
    T = Tlist(i);
    % quantize
    MWI = floor(abs(MW/T)).*sign(MW);
    % bits    
    Esep = 0;
    Jmax = log2(n)-1; Jmin = 4;
    for j = Jmax:-1:Jmin
        for q=1:3
            [selx,sely] = compute_quadsel(j,q);
            MWj = MWI(selx,sely);
            Esep = Esep + prod(size(MWj))*compute_entropy(MWj);
        end
    end
    Esep = Esep + prod(size(MWj))*compute_entropy(MWI(1:2^j,1:2^j));
    Ewhole = compute_entropy(MWI);
    nbits1(i) = Esep/n^2;
end
clf;
hh = plot([nbits(:)';nbits1(:)']',[err(:)';err(:)']'); axis('tight');
set_label('bpp','SNR');
legend('Whole','Separate');
if using_matlab()
    set(hh,'LineWidth',2);
end
%EXO
