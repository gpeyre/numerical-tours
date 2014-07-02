function PR = radon(IMG, THETA)

%% This MATLAB function takes an image matrix and vector of angles and then
%% finds the 1D projection (Radon transform) at each of the angles.  It returns
%% a matrix whose columns are the projections at each angle.
%%
%% Written by : Justin K. Romberg

% pad the image with zeros so we don't lose anything when we rotate.
[iLength, iWidth] = size(IMG);
iDiag = sqrt(iLength^2 + iWidth^2);
LengthPad = ceil(iDiag - iLength) + 2;
WidthPad = ceil(iDiag - iWidth) + 2;
padIMG = zeros(iLength+LengthPad, iWidth+WidthPad);
padIMG(ceil(LengthPad/2):(ceil(LengthPad/2)+iLength-1), ...
    ceil(WidthPad/2):(ceil(WidthPad/2)+iWidth-1)) = IMG;

% loop over the number of angles, rotate 90-theta (because we can easily sum
% if we look at stuff from the top), and then add up.  Don't perform any
% interpolation on the rotating.

n = size(padIMG,1);
x = linspace(-1,1,n);
[X1,Y1] = meshgrid(x,x);

n = length(THETA);
PR = zeros(size(padIMG,2), n);
for i = 1:n
    theta = (90-THETA(i))*pi/180;
    X = cos(theta)*X1 + -sin(theta)*Y1;
    Y = sin(theta)*X1 + cos(theta)*Y1;    
    % interpolate
    tmpimg = interp2(X1,Y1,padIMG,X,Y);
    tmpimg(isnan(tmpimg)) = 0;
    % sum
    PR(:,i) = (sum(tmpimg))';
end
