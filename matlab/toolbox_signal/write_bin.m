function write_bin(M,  name)

% write_bin - write a file to .bin format readable by scilab.
%
%   write_bin(M,  name);
% OR 
%   write_bin(name);
%
%   Convert to 8 bits graylevel images before writing.
%
%   name should not contain .bin extension (automatically added).
%
%   Copyright (c) 2008 Gabriel Peyre

if isstr(M)
    name = M;
    M = load_image(M);
    M = mean(M,3);
end

name = [name '.bin'];
fid = fopen(name, 'wb');
if fid<0
    error('Problem with file writing.');
end
M = round( rescale(M)*255 );

n = size(M,1);
p = size(M,2);
q = size(M,3);

fwrite(fid, n,'uint16');
fwrite(fid, p,'uint16');

if nb_dims(M)==3
    fwrite(fid, q,'uint16');    
end

fwrite(fid, M(:),'uint8');
fclose(fid);
