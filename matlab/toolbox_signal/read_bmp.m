function f = read_bmp(filename)

% read_bmp - read BMP files
%
%   f = read_bmp(filename);
%
%   Read uncompressed 8bits (grayscale) or 24bits (color) uncompressed BMP files.
%
%   If filename does not finish by '.bmp', then this extension is added.
%
%   Copyright (c) 2010 Gabriel Peyre

if not(strcmp(filename(end-3:end), '.bmp'))
    filename = [filename '.bmp'];
end

fid = fopen(filename, 'rb');
if fid<0
    error(['Cannot open ' filename '.']);
end
    

%%
% File Header

% file type, 2c
filetype = fread(fid, 2,'char');
filetype = char(filetype); 
% overall size, 4c 
filesize = fread(fid, 1,'int32');
% reserved
tmp = fread(fid, 1,'int32');
% offset
offs = fread(fid, 1,'int32');

%%
% Image header

headersize = fread(fid, 1,'int32');
width = abs(fread(fid, 1,'int32'));
height = abs(fread(fid, 1,'int32'));
nplanes = fread(fid, 1,'int16');    % should be 1
colordepth = fread(fid, 1,'int16');
ncol = ceil(colordepth/8);
compressiontype = fread(fid, 1,'int32'); 
if compressiontype~=0
    error('Only works for uncompressed BMP');
end
totalsize = fread(fid, 1,'int32');
resolhori = fread(fid, 1,'int32');
resolvert = fread(fid, 1,'int32');
palettesize = fread(fid, 1,'int32');
if palettesize~=0
    error('Does not works for colrmap BPM');
end
paletteimpsize = fread(fid, 1,'int32');

if ncol==1
    % BUG
    height = height+2;
end

%
linesize = ceil(ncol*width/4)*4;
f = fread(fid, linesize*height,'char');
f = reshape( f, [linesize height] );
f = f(1:ncol*width,:);
f = reshape(f, [ncol width height]);
f = permute(f, [3 2 1]);
if ncol>1
    f = f(:,:,3:-1:1);
end
f = double(f)/255;
fclose(fid);

if ncol==1
    % BUG
    f = f(3:end,:,:);
end