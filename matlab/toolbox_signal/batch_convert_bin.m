% convert image to .bin
name = 'flowers';
M = load_image(name);
write_bin(rescale(M), name);