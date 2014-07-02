function imageplot(M, str, a, b, c)

// display an image M.

if argn(2)<2
    str = [];
end

if argn(2)==5
    subplot(a,b,c);
end

M = squeeze(M);

if nb_dims(M)>2
if size(M,3)==2
    // for vector fields
    A = M; 
    M = zeros(size(M,1), size(M,2),3);
    M(:,:,1:2) = A;
end
if size(M,3)>3 
    warning('Works for less than 3 cannals');
    M = M(:,:,1:3);
end
end

M = rescale(M);
imshow(M);

if not(isempty(str))
    title(str);
end

set_axis(0);

endfunction

function imshow(img, arg2,strf)
//
// Displays image in Scilab graphic window.
//    imshow(I)      // I assumed to be 0-1
//    imshow(I,n)    // I assumed to be 0-1
//    imshow(I,[low high]) 
//    imshow(X,map)  // X assumed to be 1-N
//    imshow(RGB)    // RGB 0-1
//    imshow filename
//
// AUTHOR
//    Ricardo Fabbri  <rfabbri@if.sc.usp.br>
//    Cybernetic Vision Research Group
//    Luciano da Fontoura Costa, supervisor.
//    http://cyvision.if.sc.usp.br
//
// TODO
//
// - incorporate other Matplot parameters 
//
// $Revision: 1.18 $ $Date: 2004/02/13 20:36:42 $

[lhs, rhs] = argn(0);
if rhs == 0 then
   error('Invalid number of arguments.');
end

if ~exists('strf','local')
   strf='040'
else
   rhs=rhs-1
end

wins = 512;

// TODO
// double-buffering will be default for the next release
//
//prev_pixmap_mode = xget('pixmap');
//
//// double-buffering:
//if ~exists('sip_disable_pixbuffer') 
//   xset('pixmap',1);
//end

select type(img)
   case 1 then    // paletted or grayscale image 
      if rhs == 1 then
         [m,n]=size(img)
         xset('colormap', graycolormap(256))
         xset('wdim',n,m)
         Matplot(img*255 + 1,strf)
      else
         [m2,n2] = size(arg2)
         [m,n]=size(img)
         if n2 == 1 then
            if type(arg2) ~= 1 then    // imshow(img,ncolors)
               error("2nd argument must be a scalar.")
            end
            xset('colormap', graycolormap(arg2))
            img = img*(arg2-1)+1;
         elseif (n2==0 | n2==2) then   // [] or [n1 n2]
            if n2~=0 then
               img(img<=arg2(1))=arg2(1)
               img(img>=arg2(2))=arg2(2)
            end
            img=round(normal(img,256,1));
            xset('colormap', graycolormap(256))
         elseif (n2 == 3) then         // imshow(img,map)
            if m2 == 65536
               //
               // Scilab reserves black & white colors, so only max 65534
               // entries are permited in the colormap :-( 
               // Black and white are always accessible by m2+1 and m2+2, 
               // respectively, no matter what is the current colormap.
               //
               // I'm ashamed of this code... we must urgently improve scilab's
               // image display. Please help us if you can :-)
               //
               img2 = ind2rgb(img,arg2);
               arg2=sip_approx_true_cmap(11)  // 11 levels (11^3 colors)
               img=sip_index_true_cmap(img2,11)
            end
            xset('colormap', arg2)
         else
            error('Invalid size of 2nd argument.')
         end
         xset('wdim',wins,wins);
         Matplot(img,strf)
      end   
   case 17 then   // truecolor image
      dims=size(img)
      if dims(3) ~= 3 then
         error('RGB image must have 3rd dimension equal to 3.')
      end
      // The following works, at the cost of reduction
      // of the number of colors to 40^3 ~= 2^16
      // Scilab unfortunately doesn't work with 25bit colordepth
      // My thanks to Bruno Pincon
//      printf('This may take a while.');
      if argn(2) == 1
         nlevels = 11  // 11^3 colors only, for speed
      else
         nlevels = arg2
      end
      xset('colormap',sip_approx_true_cmap(nlevels))
      xset('wdim',wins,wins)
      Matplot(sip_index_true_cmap(img,nlevels),strf)
   case 10 then   // filename
      map=0;
      [image,map] = imread(img)
      imshow(image,map)
   else
      error('1st argument of invalid type.')
end

if xget('pixmap') == 1,  xset('wshow'), end

// Part of double-buffer default for next release:
//xset('pixmap',prev_pixmap_mode);
//

if MSDOS, xbasr(), end

endfunction


function cmap = sip_approx_true_cmap(n)
//
// There are n levels for each color channel intensity
// (each intensity being given by an integer I between 0 and n-1)
// To the "color" R,G,B (R,G,B in [0,n-1]) must correspond the
// index k= R n^2 + G n + B + 1 of the table cmap of size n^3 x 3
// and cmap(k,:) =  [R/(n-1) G/(n-1) B/(n-1)]
//
// As the max size of a cmap in scilab is 2^16-2, 
// n = 40 is the max possible (40^3 <= 2^16 - 2 < 41^3).
// 
// This function returns this colormap.
//
// ORIGINAL AUTHOR 
//	   Bruno Pincon <bruno.pincon@free.fr>
//

if argn(2)==0
   n = 40
end
nb_col = n^3
temp = (0:nb_col-1)'
cmap = zeros(nb_col,3)
q = int(temp/n^2)
cmap(:,1) = q/(n-1)
q = modulo(int(temp/n),n)
cmap(:,2) = q/(n-1)
cmap(:,3) = modulo(temp,n)/(n-1)
endfunction

function [A] = sip_index_true_cmap(Im,n)
   //
   // On input Im is a n1 x n2 x 3 hypermat describing a 
   // true color image  Im(i,j,:) giving the R-G-B of the 
   // pixel (i,j)
   //
   // On output A is a n1 x n2 matrix, A(i,j) given the 
   // index on the "true" color map of the (i,j) pixel. 
   // 
   // This new version doesn't use anymore hypermatrices
   // extraction for a gain in speed (as hypermatrix extraction 
   // is not currently too efficient in scilab); this result in
   // a less aesthetic code than before...  Also it uses
   // round in place of floor for a better color reduction
   // (from 0-255 levels to 0-39). 
   //
   // Author : Bruno Pincon
   //
   if argn(2)==1
      n = 40
   end
   dims = size(Im)
   v = round(Im("entries")*(n-1))
   m = dims(1)*dims(2)
   A = v(1:m)*n^2 + v(m+1:2*m)*n + v(2*m+1:$) + 1
   A = matrix(A,dims(1),dims(2))
endfunction


