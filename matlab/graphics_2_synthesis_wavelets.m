%% Texture Synthesis Using Wavelets
% This numerical tour explores texture synthesis using wavelets.

%%
% Image synthesis is obtained by drawing an image at random that satisfies
% some modeling constraint, that are usually learned from a given exemplar
% texture. 


perform_toolbox_installation('signal', 'general');


%% Multi-scale Texture Synthesis
% The decay of wavelet coefficients caraterize pointwise singularities in
% images and texture. Histogram equalization enable the synthesis of
% texture with singularities. This corresponds to the texture synthesis
% algorithm of Heeger and Bergen.

%%
% Load a texture.

n = 512;
name = 'texture';
M = load_image(name, n);
M = rescale( sum(M,3) );

%% 
% *For Scilab users*: you should increase the size of the memory.
% _Warning_: execute this line only once.

extend_stack_size(4);

%%
% First we compute the wavelet coefficients of the texture.
% We use a translation invariant transform.

options.ti = 1;
Jmin = 4;
MW = perform_wavelet_transf(M(:,:,1), Jmin, +1, options);

%%
% We initialize the synthesis by a random noise with the same gray values.

M1 = perform_hist_eq(randn(n,n), M);

%%
% Display.

clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Initial noise', 1,2,2);

%% 
% We also compute the wavelet transform of the noise.

MW1 = perform_wavelet_transf(M1, Jmin, +1, options);

%%
% A random texture is obtained by histogram equalization of each wavelet
% scale.

for i=1:size(MW,3)
    MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i));    
end

%%
% We retrieve the texture by inverse wavelet transform.

M1 = perform_wavelet_transf(MW1, Jmin, -1, options);

%%
% Display.

clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Initial synthesis', 1,2,2);

%EXO
%% Iterate these two steps (spatial and wavelet histogram matching) until convergence to a stable step.
% spacial matching
niter = 4;
M1 = randn(n);
clf;
for k=1:niter
    M1 = perform_hist_eq(M1, M);
    % wavelet matching
    MW1 = perform_wavelet_transf(M1, Jmin, +1, options);
    for i=1:size(MW,3)
        MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i));
    end
    M1 = perform_wavelet_transf(MW1, Jmin, -1, options);
    % display
    imageplot(M1, strcat(['Iteration ' num2str(k)]), 2,2,k);
end
%EXO

%% Multi-scale Color Texture Synthesis
% It is possible to perform color synthesis by synthesizing independantly
% each channel over a well chosen color space.

%% 
% Load a color texture.

n = 512;
M = rescale( load_image('texture', n) );

%EXO
%% Perform color texture synthesis with wavelets over the RGB space.
niter = 3;
M1 = randn(n,n,3);
for c=1:3
    MW = perform_wavelet_transf(M(:,:,c), Jmin, +1, options);
    for k=1:niter
        M1(:,:,c) = perform_hist_eq(M1(:,:,c), M(:,:,c));
        % wavelet matching
        MW1 = perform_wavelet_transf(M1(:,:,c), Jmin, +1, options);
        for i=1:size(MW,3)
            MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i));
        end
        M1(:,:,c) = perform_wavelet_transf(MW1, Jmin, -1, options);
    end
end
% Display.
clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Synthesized', 1,2,2);
%EXO

%EXO
%% Try with other color spaces, for instance PCA adapte space.
%EXO

%% Multi-dimensional Color Equalization
% To maintain color consistency, it is possible to use a color
% equalization.

%%
% Initial image.

M1 = randn(n,n,3);

%% 
% A simple (but not very acurate) method to perform in performing
% independant channel equalization over randomized color space. 
% This needs to be repeated several time to converge to a real matching.

%%
% Compute a random 3x3 orthogonal matrix.

[U,R] = qr(randn(3));

%%
% Perform the change of color space.

d = reshape(M,[n^2 3])*U;
d1 = reshape(M1,[n^2 3])*U;

%%
% Perform the equalization

for c=1:3
    d1(:,c) = perform_hist_eq(d1(:,c),d(:,c));
end

%%
% Perform the inverse change of color space.

M1 = reshape(d1*U',[n n 3]);

%%
% Compares the histogram of the R channel. You can see that the match is
% imperfect.

m = M(:,:,1); m1 = M1(:,:,1);
clf;
subplot(2,1,1);
hist(m(:),50); title('Original');
subplot(2,1,2);
hist(clamp(m1(:)),50); title('Matched'); 

M1 = randn(n,n,3);
for i=1:3
    M1(:,:,i) = perform_hist_eq(M1(:,:,i), M(:,:,i));
end

%EXO
%% Perform iteratively the randomized matching. Plot the decay of the
%% mathing error.
M1 = randn(n,n,3);
niter1 = 200;
delta = [];
for k=1:niter1
    [U,R] = qr(randn(3));
    d = reshape(M,[n^2 3])*U;
    d1 = reshape(M1,[n^2 3])*U;
    delta(k) = 0;
    for c=1:3
        [tmp,I] = sort(d(:,c)); 
        [tmp,I1] = sort(d1(:,c)); 
        delta(k) = delta(k) + norm(d1(I1,c) - d(I,c))^2;
        d1(I1,c) = d(I,c);
        % d1(:,c) = perform_hist_eq(d1(:,c),d(:,c));
    end
    M1old = M1;
    M1 = reshape(d1*U',[n n 3]);
    % delta(k) = norm(M1(:)-M1old(:));
end
clf;
plot(log10(delta/norm(M1(:))), '.-');
axis('tight');
%EXO

%%
% Display the histograms of the R channels. The match is not perfect, 
% but it is better than with a single projection.

m = M(:,:,1); m1 = M1(:,:,1);
clf;
subplot(2,1,1);
hist(m(:),50); title('Original');
subplot(2,1,2);
hist(clamp(m1(:)),50); title('Matched');

%%
% Display the equalized color image;

clf;
imageplot(M, 'Image', 1,2,1);
imageplot(M1, 'Equalized', 1,2,2);

%EXO
%% Perform color texture synthesis with wavelets using this color histogram
%% matching at each iteration.
niter = 5;
niter1 = 5;
M1 = randn(n,n,3);
M1 = perform_color_matching(M1,M,niter1);
% precompute the wavelet transform of the image
for c=1:3
    MW(:,:,:,c) = perform_wavelet_transf(M(:,:,c), Jmin, +1, options);
end
Msvg1 = {}; Msvg = {};
for k=1:niter
    % wavelet matching
    for c=1:3
        % transform
        MW1 = perform_wavelet_transf(M1(:,:,c), Jmin, +1, options);
        % equalize
        for i=1:size(MW1,3)
            MW1(:,:,i) = perform_hist_eq(MW1(:,:,i), MW(:,:,i,c));
        end
        M1(:,:,c) = perform_wavelet_transf(MW1, Jmin, -1, options);
    end
    % spacial matching
    M1 = perform_color_matching(M1,M,niter1);
    Msvg{end+1} = M1;
end
% Display.
clf;
imageplot(M, 'Exemplar', 1,2,1);
imageplot(M1, 'Synthesized', 1,2,2);
%EXO

