%% Stationary Dynamic Texture Synthesis
% This tour explores the analysis and synthesis of stationary dynamic
% textures using Gaussian models.

perform_toolbox_installation('signal', 'general');

%% Video Loading
% A color video is a 4-D array of values of size \((n_1,n_2,p,3)\) where
% \(N=n_1 \times n_2\) is the number of pixels in the video and \(p\) of frames.

%%
% Load a video from a |.gif| file, which requires conversion from indexed color to RGB colors. 

name = 'smoke';
name = 'fire';
[X, map] = imread([name '.gif'], 'frames', 'all');
f = [];
for i=1:size(X,4)
    f(:,:,:,i) = ind2rgb(X(:,:,1,i),map);
end

%%
% Modify the video (here time reverse it).

g = f(:,:,:,end:-1:1);

%%
% Save the video as a |.gif| file.
% Compute the quantized color map from the first frame.

X = [];
[X(:,:,1,1),map] = rgb2ind(g(:,:,:,1),128*2);
for i=2:size(g,4)
    [X(:,:,1,i),map] = rgb2ind(g(:,:,:,i),map);    
end
imwrite(X+1, map, ['../html/graphics_9_dyntextures/video1.gif'], 'DelayTime',1/10, 'loopcount', Inf);

%%
% <html>
% <div align="center">
% <img src="video1.gif" width="256"/>
% </div>
% </html>