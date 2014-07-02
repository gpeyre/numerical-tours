%% Texture Synthesis with PDEs
% This numerical tour explores image synthesis using diffusion equation.
% You might first have a look at the numerical tours on Heat diffusion and 
% Total variation minimization that introduces these two notions.

perform_toolbox_installation('signal', 'general');


%% Heat Flow Synthesis
% It is possible to synthesis a cloudy image by a simple smoothing of a
% noise. Histogram equalization helps to maintain the contrast.

n = 256;

%%
% We blur a noisy image to have a cloud texture.

M1 = perform_blurring(randn(n,n), 15);

%%
% We impose a flat histogram to enhance the contrast.

x = linspace(0,1, n*n);
M2 = perform_hist_eq(M1,x);
% display
clf;
imageplot(M1, 'Original', 1,2,1);
imageplot(M2, 'Equalized', 1,2,2);

%EXO
%% Perform a synthesis by running a heat diffusion, starting with a random
%% noise a time |T=0|. At each step of the diffusion, perform an histogram
%% equation to keep the contrast of the texture.
T = 20;
% step size
tau = .2;
% number of iteration
niter = round(T/tau);
Mheat = perform_hist_eq(randn(n), x);
clf; k = 0;
for i=1:niter
    % compute the gradient
    G = div(grad(Mheat));
    % descent 
    Mheat = perform_hist_eq(Mheat + tau*G, x);
    if mod(i,floor(niter/4))==0
        k = k+1;
        imageplot(Mheat, strcat(['T=' num2str(T*k/4,3)]), 2,2,k );
    end
end
%EXO

%% Total Variation Constraints
% The total variation minimization reduce the number of edges and is
% usually used as a denoising method. It can also be used to perform
% synthesis of singular images.

%%
% The total variation roughly measures the number of edges present in an
% image. 

M = rescale( load_image('lena', n) );
Gr = grad(M);
tv = sum(sum( sqrt(sum(Gr.^2,3)), 2 ), 1);
disp(strcat(['Total variation of M = ' num2str(tv) '.']));

%EXO
%% Starting from an initial noise image, perform a total variation
%% minimization. At each step of the descent, perform an histogram
%% equalization so that the texture has a flat histogram.
T = 3;
% step size
tau = .005;
% avoid instabilities
epsilon = 1e-5;
% number of iteration
niter = round(T/tau);
Mtv = perform_hist_eq(randn(n), x);
clf; k = 0;
for i=1:niter
    % compute the gradient
    G = grad(Mtv);
    d = max(epsilon, sqrt(sum(G.^2,3)));
    G = div( G ./ repmat(d, [1 1 2]) );
    % descent 
    Mtv = perform_hist_eq( Mtv + tau*G, M);
    if mod(i,floor(niter/4))==0
        k = k+1;
        Gr = grad(Mtv);
        tv = sum(sum( sqrt(sum(Gr.^2,3)), 2 ), 1);
        imageplot(Mtv, strcat(['T=' num2str(T*k/4,3) ', TV=' num2str(tv)]), 2,2,k );
    end
end
%EXO

%EXO
%% Perfrom a synthesis that mixes both TV minimization (to reduce the TV
%% norm)
%% and wavelet histogram equalization (to control the distribution of singularities). Stop the iterations when the
%% synthesized image has the same TV norm as the original one.
%EXO


