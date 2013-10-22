function write_gim(M, filename, save_png)

% write_gim - write a SGIM in .gim file format (i.e. single precision binary file).
%
%   write_gim(M, filename, save_png);
%
%   If filename=='mymesh' then it will produce files 'filename.gim' and 'filename.png'.
%   If save_png==1, it will save as an image file quantified on 8 bits per channels.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    save_png = 1;
end

% save binary file
str = [filename, '.gim'];
fid = fopen(str,'wb');
fwrite(fid,M,'float');
fclose(fid);

% be sure that image lie in [0,1], so clamp
if mmax(M)>1.1 || mmin(M)<-0.1
    warning('Geometry images should lie in [0,1]');
end
% M = rescale(M);
M = max(M, 0);
M = min(M, 1);


% save PNG
if save_png
    str = [filename, '.png'];
    warning off;
    imwrite( M, str, 'png' );
    warning on;
end