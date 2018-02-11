function write_bmp(M, name)

// write_bmp - write a 8 bit bw binary file
//
//   M = write_bmp(name);
//
//   Important: M should be in [0,1].
//   The size of M should be a multiple of 4.
//
//   Copyright (c) 2008 Gabriel Peyre

if isempty((strfind(name, '.')))
    name = [name '.bmp'];
end

fid = mopen(name, 'wb');
if fid<0
    error(['File ' name ' does not exists.']);
end

M = M(size(M,1):-1:1,:);
M=M';
M = round( clamp(M*255, 0,255) );

toReachMOf4 = write_bmp_header(fid, size(M,1), size(M,2));

// write data
M = repmat(M, [1 1 3]);
M = permute(M, [3 1 2]);
mtlb_fwrite(fid,M(:),'us');  % ubit8 -> us

//write filler bites to reach multiple of 4
extra = zeros(1,toReachMOf4);
mtlb_fwrite(fid,extra,'us');

mclose(fid);

endfunction

/////
function toReachMOf4 = write_bmp_header(fid ,width, height)

bitsPerPixel = 24;
compression = 0;

bytesPerRow = width*bitsPerPixel/8;
//Number of bytes written to each row must be a mutliple of 4. Append 0's
// to end to make it fit
if mod(bytesPerRow,4) == 0
	toReachMOf4 = 0;
else
	toReachMOf4 = 4 - mod(bytesPerRow,4);
end

sizeOfFileInfo = 14;
sizeOfInfoHeader = 40;
sizeOfRGBQuadArray = 0;

offsetOfImage = sizeOfFileInfo + sizeOfInfoHeader + sizeOfRGBQuadArray;
sizeOfFile = offsetOfImage + bitsPerPixel/8*width*height;

// fid = fopen(a_name,'w');

//writes each char of the string as a ubit8
// bytes 0x 0-1
mtlb_fwrite(fid,'BM','us');

//file size
// bytes 0x 2-5
mtlb_fwrite(fid,sizeOfFile,'ubit32');

// set to zero
// bytes 0x 6-9
mtlb_fwrite(fid,[0],'ubit32');

// offset of bitmap data
// bytes 0x A-D
mtlb_fwrite(fid,offsetOfImage,'ubit32');

// size of info header in bytes
// b 0x E-11
mtlb_fwrite(fid,sizeOfInfoHeader,'ubit32');

//width of bitmap in pixels
// b 0x 12-15
mtlb_fwrite(fid,width,'ubit32');

//height of bitmap in pixels
// b 0x 16-19
mtlb_fwrite(fid,height,'ubit32');

// set to one
// b 0x 1A-1B
mtlb_fwrite(fid,1,'ubit16');

//bits per pixel
// b 0x 1C-1D
mtlb_fwrite(fid,bitsPerPixel,'ubit16');

// file compression
// b 0x 1E-21
mtlb_fwrite(fid,compression,'ui');

//size of image data or, if no compression, 0
// b 0x 22-25
mtlb_fwrite(fid,0,'ui');

// bits/meter or zero (wide, high)
// b 0x 26-29,2A-2D
mtlb_fwrite(fid,[0,0],'ui');

//	# of colors used or 0
// b 0x 2E-31
mtlb_fwrite(fid,0,'ui');

// numbers of colors that are important or zero if all are
// b 0x 31-34
mtlb_fwrite(fid,0,'ui');

//fill the space between now and the offset
//Somehwere along the lines, I lost count. Previosuly, this was set to 0x35,
//instead of 0x36, which was leading to every bite being offset by 1, which
//caused problems with colormapping, and lead to a broder around the image.
//At some point I'll find out where I lost count
//TODO find out where I lost count.
diff = offsetOfImage - hex2dec('36');
x = rand(1,diff);
mtlb_fwrite(fid,x,'us');


endfunction