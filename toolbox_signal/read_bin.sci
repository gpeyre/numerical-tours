function M = read_bin(name,options)

// read_bin - read bin file
//
//   M = read_bin(name,options);
//
//   Set options.nbdims = 3 for 3D data.
//
//   Copyright (c) 2008 Gabriel Peyre

options.null = 0;

if isempty(strfind(name, '.'))
    name = [name '.bin'];
end
image_rep = getoptions(options, 'image_rep', './toolbox_signal/');
nbdims = getoptions(options, 'nbdims', 2);

name = strcat([image_rep name]);

fid = mopen(name, 'rb');

if fid<0
    error('Problem with file loading.');
end


n = mget(1, 'us', fid);
p = mget(1, 'us', fid);
q = 1;
if nbdims==3
    q = mget(1, 'us', fid);
end
M = mget(n*p*q, 'uc', fid);
mclose(fid);

if length(M(:))~=n*p*q
    error('Problem with file loading, wrong size.');
end

if q==1
    M = matrix(double(M(:)), n,p );
else
    M = matrix(double(M(:)), n,p,q );
end

M = squeeze(M);

endfunction
