%% Mathematical Morphology
% This numerical tour explores mathematical morphology of binary images.

perform_toolbox_installation('signal', 'general');

%% Binary Images and Structuring Element
% Here we process binary images using local operator defined using a
% structuring element, which is here chosen to be a discrete disk of
% varying radius.

%%
% Load an image

n = 256;
M = rescale( load_image('cortex',n) );

%%
% Display.

clf;
imageplot(M);

%%
% Make it binary.

M = double(M>.45);

%%
% Display.

clf;
imageplot(M);

%%
% Round structuring element.

wmax = 7;
[Y,X] = meshgrid(-wmax:wmax, -wmax:wmax);
normalize = @(x)x/sum(x(:));
strel = @(w)normalize( double( X.^2+Y.^2<=w^2 ) );

%EXO
%% Display structuring elements of increasing sizes.
clf;i = 1;
for w=[1 1.5 3 7]
    imageplot(strel(w), ['w=' num2str(w)], 1,4, i);
    i = i+1;
end
%EXO

%% Dillation
% A dilation corresponds to take the maximum value of the image aroung each pixel, 
% in a region equal to the structuring element.

%%
% It can be implemented using a convolution with the structuring element
% followed by a thresholding.

dillation=@(x,w)double(perform_convolution(x,strel(w))>0);
Md = dillation(M,2);

%%
% Display.

clf;
imageplot(Md);


%EXO
%% Test with structing elements of increasing size.
clf;
i = 0;
for w = [1 1.5 2 4]
    i = i+1;
    imageplot(dillation(M,w), ['w=' num2str(w)], 2,2,i);
end
%EXO


%% Errosion
% An errosion corresponds to take the maximum value of the image aroung each pixel, 
% in a region equal to the structuring element.

%%
% It can be implemented using a convolution with the structuring element
% followed by a thresholding.

errosion=@(x,w)double( perform_convolution(x,strel(w))>=.999 );
Me = errosion(M,2);


%%
% Display.

clf;
imageplot(Me);


%EXO
%% Test with structing elements of increasing size.
clf;
i = 0;
for w = [1 1.5 2 4]
    i = i+1;
    imageplot(errosion(M,w), ['w=' num2str(w)], 2,2,i);
end
%EXO


%% Opening
% An opening smooth the boundary of object (and remove small object) by
% performing an errosion and then a dillation.

%%
% Define a shortcut.

opening = @(x,w)dillation(errosion(x,w),w);

%%
% Perform the opening, here using a very small disk.

w = 1;
Mo = opening(M,w);

%%
% Display.

clf;
imageplot(Mo);

%EXO
%% Test with structing elements of increasing size.
clf;
i = 0;
for w = [1 1.5 2 4]
    i = i+1;
    imageplot(opening(M,w), ['w=' num2str(w)], 2,2,i);
end
%EXO
