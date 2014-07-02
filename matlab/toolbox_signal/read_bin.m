function M = read_bin(name,options)

% read_bin - read bin file
%
%   M = read_bin(name,options);
%
%   Set options.nbdims = 3 for 3D data.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;

if isempty(strfind(name, '.'))
    name = [name '.bin'];
end
image_rep = getoptions(options, 'image_rep', './toolbox_signal/');


nbdims = getoptions(options, 'nbdims', 2);

name = [image_rep name];

fid = fopen(name, 'rb');

if fid<0
    error('Problem with file loading.');
end


n = fread(fid, 1,'uint16');
p = fread(fid, 1,'uint16');
q = 1;
if nbdims==3
    q = fread(fid, 1,'uint16');
end

M = fread(fid, n*p*q, 'uint8');
fclose(fid);

if length(M)~=n*p*q
    error('Problem with file loading.');
end
M = reshape(M, [n, p, q]);
