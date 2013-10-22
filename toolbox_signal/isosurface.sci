function isosurface(F,t)

// isosurface - draw an isosurface
//
//  isosurface(F,t);
//
//  Uses the code 
//  Color graphics, 3d rendering and data manipulation routines for Scilab 
//  by Enrico Segre
//  http://www.weizmann.ac.il/home/fesegre/scistuff.html


f = gcf()
f.color_map = graycolormap(1024);
drawlater();

x = 1:size(F,1);
y = 1:size(F,2);
z = 1:size(F,3);

[xx,yy,zz]=isosurf3d(x,y,z,F,t);
shadesurf(xx,yy,zz);
drawnow();

endfunction;

function [xx,yy,zz]=isosurf3d(x,y,z,s,s0)
// 
//   xx,yy,zz= (nf,3)  facelet vertices
//   s= real hypermat([length(x),length(y),length(z)])
//   s0 scalar
//
//   hexaedral FEM syntax:
//   [xx,yy,zz]=isosurf3d(dcorvg,kvert,field,f0)
//
//     dcorvg=(nvert,3) coordinates of the vertices of the cells
//     kvert=(ncells,8) list of the vertices of each cell (in a proper order)
//     field=(nvert,1) values of the scalar field at the vertex points
//     f0 scalar 
//
// example: see the demo below.
//
// Plot with:
//    xbasc(); plot3d1({xx;xx(1,:)},{yy;yy(1,:)},{zz;zz(1,:)})
// or
//    xbasc(); shadesurf(xx,yy,zz)
//
// note: this implementation is rather slow, surely due to the
//       nested loops and hypermat operations, which I don't know
//       how to vectorize more readably with scilab syntax (in
//       principle the algorithm is completely parallelizable). 
//       Timings are data dependent, as the work depends on
//       the number of facelets found.

//t=0

[lhs,rhs]=argn(0);
if rhs==0 then
//demo
    disp "demo of isosurf3d(x,y,z,s,s0)"
    disp " "
    democomm=[
    "  nx=20; ny=20; nz=20; s=hypermat([nx,ny,nz]);"
    "  x=linspace(-4.5,4.5,nx); xq=x''*ones(1,ny);"
    "  y=linspace(-4.5,4.5,ny); yq=ones(nx,1)*y;"
    "  z=linspace(-4.5,4.5,nz);"
    "  for i=1:nz; s(:,:,i)=1../((xq-2).^2+yq.^2+(z(i)-2)^2)+"+..
    "1../((xq+2).^2+yq.^2+(z(i)+2)^2); end"
    "  timer();[xx,yy,zz]=isosurf3d(x,y,z,-s,-.18); t=timer();"
    ]
    write(%io(2),democomm); execstr(democomm)
// the -0.2 surface gives a good example of the saddle point hole, btw
    aa=string(size(xx,2))+" triangles, "+string(t)+" sec --> "+..
       string(size(xx,2)/t)+" triang./sec"
    disp "  xbasc(); shadesurf(xx,yy,zz,1,0,60,60,''x@y@z'',[1 6 4])"
    xbasc();  shadesurf(xx,yy,zz,1,0,60,60,'x@y@z',[1 6 4])
    xx=aa;
    return
end;

if rhs==2 & type(x)==17 then
   s=x; s0=y; nx=size(s,1); ny=size(s,2); nz=size(s,3); 
   x=1:nx; y=1:ny; z=1:nz;
end 

if exists('s','local') & type(s)==17 then 
   datastruct='h';
   nx=size(s,1); ny=size(s,2); nz=size(s,3); 
   x=matrix(x,1,size(x,'*')); y=matrix(y,1,size(y,'*')); 
   z=matrix(z,1,size(z,'*'));
   if length(x)<nx | length(y)<ny | length(z)<nz then
     disp 'check the dimensions of the arguments!'; return
   end 
else
   datastruct='c'; 
   dcorvg=x; kvert=y; field=z; f0=s; 
   nvert=size(dcorvg,1); ncell=size(kvert,1)
   if size(field)<>[nvert,1] then
      disp 'wrong dimensions of the scalar field!'; return
   end
end


// Now, we have to generate a set of facelets which approximate
// the isosurface sought. My approach is to scan the hypermatrix
// cell by cell, and a) to identify the pattern of vertices > and
// < s0 (there are 23 possible patterns, which can appear in one
// of 48 possible orientations), then b) to generate a set of
// triangular facelets which represent the isosurface inside the cell.
// Ideally, the function inside the cell would be better
// approximated by shape functions. Shape functions evaluating
// to 1 in one vertex and to 0 in the other 7 would be of higher
// than linear order, and their isosurfacelets wouldn't be
// plane. My (to some extent arbitrary) approximation as
// triangular facelets whose vertices are the linear solutions
// along the sides of the cube, should at least be easier.
// I've found this method to work almost always, the only 
// exception being that it leaves a hole in some 2-cell saddle
// points (of type [[0 1;1 0],[1 0;0 1]]). Treating that would
// require to consider couples of adjacent cells, which is more
// complicate. Alternative approaches would of course be possible
// (e.g., only rectangles parallel to the coordinate planes
//  on the dual lattice,...).
//
//         8----7
//        /|   /|
//       5----6 |     reference order of vertices in the base
//       | 4--|-3     cell
//       |/   |/
//       1----2
//
//the following is obtained with [pp,par]=cubeperm() below, 
// and snippets like
//
//for i=1:48; write(%io(2),string(pp(i,1))+' '+string(pp(i,2))+' '+..
//          string(pp(i,3))+' '+string(pp(i,4))+' '+..
//          string(pp(i,5))+' '+string(pp(i,6))+' '+..
//          string(pp(i,7))+' '+string(pp(i,8))+'; ..' ); end
// time savings are minimal, but one never knows if they don't turn
// out useful
 
 
pp=[8 7 6 5 4 3 2 1; 8 7 3 4 5 6 2 1; 8 5 6 7 4 1 2 3; ..
    8 5 1 4 7 6 2 3; 8 4 3 7 5 1 2 6; 8 4 1 5 7 3 2 6; ..
    7 8 5 6 3 4 1 2; 7 8 4 3 6 5 1 2; 7 6 5 8 3 2 1 4; ..
    7 6 2 3 8 5 1 4; 7 3 4 8 6 2 1 5; 7 3 2 6 8 4 1 5; ..
    6 7 8 5 2 3 4 1; 6 7 3 2 5 8 4 1; 6 5 8 7 2 1 4 3; ..
    6 5 1 2 7 8 4 3; 6 2 3 7 5 1 4 8; 6 2 1 5 7 3 4 8; ..
    5 8 7 6 1 4 3 2; 5 8 4 1 6 7 3 2; 5 6 7 8 1 2 3 4; ..
    5 6 2 1 8 7 3 4; 5 1 4 8 6 2 3 7; 5 1 2 6 8 4 3 7; ..
    4 8 7 3 1 5 6 2; 4 8 5 1 3 7 6 2; 4 3 7 8 1 2 6 5; ..
    4 3 2 1 8 7 6 5; 4 1 5 8 3 2 6 7; 4 1 2 3 8 5 6 7; ..
    3 7 8 4 2 6 5 1; 3 7 6 2 4 8 5 1; 3 4 8 7 2 1 5 6; ..
    3 4 1 2 7 8 5 6; 3 2 6 7 4 1 5 8; 3 2 1 4 7 6 5 8; ..
    2 6 7 3 1 5 8 4; 2 6 5 1 3 7 8 4; 2 3 7 6 1 4 8 5; ..
    2 3 4 1 6 7 8 5; 2 1 5 6 3 4 8 7; 2 1 4 3 6 5 8 7; ..
    1 5 8 4 2 6 7 3; 1 5 6 2 4 8 7 3; 1 4 8 5 2 3 7 6; ..
    1 4 3 2 5 8 7 6; 1 2 6 5 4 3 7 8; 1 2 3 4 5 6 7 8];

ppar=[1 -1 -1  1  1 -1 -1  1  1 -1 -1  1 -1  1  1 -1 -1  1 ..
      1 -1 -1  1  1 -1 -1  1  1 -1 -1  1  1 -1 -1 ..
      1  1 -1  1 -1 -1  1  1 -1 -1  1  1 -1 -1  1];

// "parity" of the permutation of 1:8 above. This is important
// to discriminate the inside and the outside of the isosurface element

[ss,ppinv]=sort(-pp,'c'); 
// this takes up a negligible time, it's not worth to write
//  also ppinv as data


//// base patterns of 0 and 1 on the vertices
//b=zeros(23,8); nft=zeros(23);
//b(1,:)= [0 0 0 0 0 0 0 0];  nft(1)=0;    // 0 t
//b(2,:)= [1 0 0 0 0 0 0 0];  nft(2)=1;    // 1 t
//b(3,:)= [1 1 0 0 0 0 0 0];  nft(3)=2;    // 2 t
//b(4,:)= [1 0 1 0 0 0 0 0];  nft(4)=2;    // 2 t
//b(5,:)= [1 1 1 0 0 0 0 0];  nft(5)=3;    // 3 t
//b(6,:)= [1 0 0 0 0 0 1 0];  nft(6)=2;    // 2 t
//b(7,:)= [1 1 0 0 0 0 1 0];  nft(7)=3;    // 3 t
//b(8,:)= [1 0 1 0 0 0 0 1];  nft(8)=3;    // 3 t
//b(9,:)= [1 1 1 1 0 0 0 0];  nft(9)=2;    // 2 t
//b(10,:)=[1 1 0 0 0 0 1 1];  nft(10)=4;    // 4 t
//b(11,:)=[1 0 1 0 1 0 1 0];  nft(11)=4;    // 4 t
//b(12,:)=[1 0 1 0 0 1 0 1];  nft(12)=4;    // 4 t
//b(13,:)=[1 1 1 0 1 0 0 0];  nft(13)=4;    // 4 t
//b(14,:)=[1 1 1 0 0 1 0 0];  nft(14)=4;    // 4 t 
//b(15,:)=[1 1 1 0 0 0 0 1];  nft(15)=4;    // 4 t
//b(16:23,:)=abs(1-b(8:-1:1,:)); nft(16:23)=nft(8:-1:1);
//bb=b==1; 
//// now that patterns are inlined as data, the array b isn't
//// be needed anymore

nft=[0 1 2 2 3 2 3 3 2 4 4 4 4 4 4 3 3 2 3 2 2 1 0]
// number of triangles for each of the 23 basic patterns

//bp=zeros(23*48,8)==1;
//for l=1:23; bp(l:23:$,:)=matrix(bb(l,pp),48,8); end;
//bp1=bp(:,1); bp2=bp(:,2); bp3=bp(:,3); bp4=bp(:,4); 
//bp5=bp(:,5); bp6=bp(:,6); bp7=bp(:,7); bp8=bp(:,8); 
//// these bp1..bp8 make up the columns of the matrix of
//// all patterns permuted in all possible ways


// definition of the triangular facelets for each pattern:
//  (all but full cell/void cell, which are skipped)
// The entries of vtrian identify, for each triangle, the 
//  three sides of the cell on which a vertex lies. Each of 
//  these sides is identified by the couple of vertices of the
//  cell among which it runs. When the pattern requires less
//  than 4 triangles, the extra entries of vtrian are ones.
// vtrian=ones(45,6*4);
// for l=-23:23
//        nt=0
//        select abs(l)
//          case 2 then 
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,4,1,5];  nt=nt+1 ,
//          case 3 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,3,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,6,2,3];  nt=nt+1 ,
//          case 4 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,4,1,5];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,7,3,4,2,3];  nt=nt+1 ,
//          case 5 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,1,5,2,6];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,2,6,3,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,7,3,4,2,6];  nt=nt+1 ,
//          case 6 then
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,4,1,5];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [6,7,7,8,3,7];  nt=nt+1 ,
//          case 7 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,3,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,6,2,3];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [6,7,7,8,3,7];  nt=nt+1 ,
//          case 8 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,4,1,5];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,7,3,4,2,3];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [5,8,7,8,4,8];  nt=nt+1 ,
//          case 9 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,6,3,7];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,3,7,4,8];  nt=nt+1 ,
//          case 10 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,3,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,6,2,3];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [4,8,6,7,5,8];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [6,7,4,8,3,7];  nt=nt+1 ,
//          case 11 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,4,5,6];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,5,8,5,6];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,3,7,8,3,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,3,6,7,7,8];  nt=nt+1 ,
//          case 12 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,4,1,5];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,7,3,4,2,3];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [4,8,7,8,5,8];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,6,5,6,6,7];  nt=nt+1 ,
//          case 13 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,5,8,5,6];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [5,6,2,6,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,6,3,4,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,6,3,7,3,4];  nt=nt+1 ,
//          case 14 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,5,6,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [5,6,3,4,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [5,6,6,7,3,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,4,6,7,3,7];  nt=nt+1 ,
//          case 15 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,1,5,2,6];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,2,6,3,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [4,8,7,8,5,8];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,6,3,7,3,4];  nt=nt+1 ,
//          case 16 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,5,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,4,3,7,2,3];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [4,8,7,8,5,8];  nt=nt+1 ,
//          case 17 then,
//             vtrian(l+23,6*nt+(1:6)) = [2,3,1,5,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,5,2,6,2,3];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,7,7,8,6,7];  nt=nt+1 ,
//          case 18 then
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,5,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,7,7,8,6,7];  nt=nt+1 ,
//          case 19 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,2,6,1,5];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [1,4,3,4,2,6];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,4,3,7,2,6];  nt=nt+1 ,
//          case 20 then,
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,5,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [3,4,3,7,2,3];  nt=nt+1 ,
//          case 21 then,
//             vtrian(l+23,6*nt+(1:6)) = [2,3,1,5,1,4];  nt=nt+1 ,
//             vtrian(l+23,6*nt+(1:6)) = [2,3,2,6,1,5];  nt=nt+1 ,
//          case 22 then 
//             vtrian(l+23,6*nt+(1:6)) = [1,2,1,5,1,4];  nt=nt+1 ,
//        end
// // switch the order of vertices if the parity of the transformation was
// // negative - this way the 'upper' of any facelet points consistently
// //  in the direction where the function increases
//        if  l<0 & l>-23 then
//           aa=vtrian(l+23,[1,2,7,8,13,14,19,20]); 
//           vtrian(l+23,[1,2,7,8,13,14,19,20])=...
//             vtrian(l+23,[5,6,11,12,17,18,23,24]);
//           vtrian(l+23,[5,6,11,12,17,18,23,24])=aa
//        end
// end
// to produce the output for the inlining
//   for i=1:45;printf('%i ',vtrian(i,:)');end

vtrian = [ 
1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 2 3 1 5 2 6 2 3 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 1 2 2 3 3 7 3 4 1 4 1 5 1 2 1 4 1 5 1 2;
1 5 2 6 1 4 2 6 3 4 1 4 2 6 3 7 3 4 1 4 1 5 1 2;
1 4 1 5 1 2 6 7 7 8 3 7 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 2 3 2 3 2 6 1 5 6 7 7 8 3 7 1 4 1 5 1 2;
1 4 1 5 1 2 2 3 3 7 3 4 5 8 7 8 4 8 1 4 1 5 1 2;
2 6 1 5 1 4 3 4 2 6 1 4 5 8 7 8 4 8 3 4 3 7 2 6;
1 4 5 6 1 5 1 4 3 4 5 6 3 4 6 7 5 6 3 7 6 7 3 4;
5 6 5 8 1 4 1 4 2 6 5 6 1 4 3 4 2 6 3 4 3 7 2 6;
1 5 1 4 1 2 2 3 3 4 3 7 5 8 7 8 4 8 6 7 5 6 2 6;
5 6 1 4 1 2 5 6 5 8 1 4 3 4 7 8 2 3 7 8 6 7 2 3;
1 4 2 3 1 5 2 3 2 6 1 5 5 8 6 7 4 8 3 7 4 8 6 7;
3 7 2 6 1 5 4 8 3 7 1 5 1 4 1 5 1 2 1 4 1 5 1 2;
1 5 1 4 1 2 2 3 3 4 3 7 4 8 7 8 5 8 1 4 1 5 1 2;
1 4 2 3 1 5 2 3 2 6 1 5 3 7 7 8 6 7 1 4 1 5 1 2;
1 5 1 4 1 2 3 7 7 8 6 7 1 4 1 5 1 2 1 4 1 5 1 2;
2 6 1 5 1 4 3 4 2 6 1 4 2 6 3 4 3 7 1 4 1 5 1 2;
1 5 1 4 1 2 2 3 3 4 3 7 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 2 3 1 5 2 3 2 6 1 5 1 4 1 5 1 2 1 4 1 5 1 2;
1 5 1 4 1 2 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2;
1 2 1 4 1 5 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2;
1 5 2 3 1 4 1 5 2 6 2 3 1 4 1 5 1 2 1 4 1 5 1 2;
1 2 1 4 1 5 3 7 3 4 2 3 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 1 5 2 6 1 4 2 6 3 4 3 7 3 4 2 6 1 4 1 5 1 2;
1 2 1 4 1 5 6 7 7 8 3 7 1 4 1 5 1 2 1 4 1 5 1 2;
1 5 2 3 1 4 1 5 2 6 2 3 6 7 7 8 3 7 1 4 1 5 1 2;
1 2 1 4 1 5 3 7 3 4 2 3 5 8 7 8 4 8 1 4 1 5 1 2;
1 5 2 6 3 7 1 5 3 7 4 8 1 4 1 5 1 2 1 4 1 5 1 2;
1 5 2 3 1 4 1 5 2 6 2 3 4 8 6 7 5 8 6 7 4 8 3 7;
1 2 1 4 5 6 1 4 5 8 5 6 2 3 7 8 3 4 2 3 6 7 7 8;
1 2 1 4 1 5 3 7 3 4 2 3 4 8 7 8 5 8 2 6 5 6 6 7;
1 4 5 8 5 6 5 6 2 6 1 4 2 6 3 4 1 4 2 6 3 7 3 4;
1 5 5 6 1 4 5 6 3 4 1 4 5 6 6 7 3 4 3 4 6 7 3 7;
1 4 1 5 2 6 1 4 2 6 3 4 4 8 7 8 5 8 2 6 3 7 3 4;
1 2 1 5 1 4 3 4 3 7 2 3 4 8 7 8 5 8 1 4 1 5 1 2;
2 3 1 5 1 4 1 5 2 6 2 3 3 7 7 8 6 7 1 4 1 5 1 2;
1 2 1 5 1 4 3 7 7 8 6 7 1 4 1 5 1 2 1 4 1 5 1 2;
1 4 2 6 1 5 1 4 3 4 2 6 3 4 3 7 2 6 1 4 1 5 1 2;
1 2 1 5 1 4 3 4 3 7 2 3 1 4 1 5 1 2 1 4 1 5 1 2;
2 3 1 5 1 4 2 3 2 6 1 5 1 4 1 5 1 2 1 4 1 5 1 2;
1 2 1 5 1 4 1 4 1 5 1 2 1 4 1 5 1 2 1 4 1 5 1 2 ] 
// here I've set the extra vertices beyond 3*nt as repeated 
// 1 4 1 5 1 2 to avoid divisions by zero in vectorized 
// interpolation (doesn't really work yet)

////this is the snippet for generating the magic numbers below.
/// Performance improves if they are inlined as data.
//pfa=zeros(1,256); bfa=pfa; 

pow2=[128 64 32 16 8 4 2 1]';

//for l=0:255 
//       d2b=zeros(1,8)
//       for i=1:8
//         d2b(i)=int((l-d2b*pow2)/pow2(i))
//       end
//       sb8=d2b==1
//       pbf=find(bp1==sb8(1) & bp2==sb8(2) & bp3==sb8(3)..
//              & bp4==sb8(4) & bp5==sb8(5) & bp6==sb8(6)..
//              & bp7==sb8(7) & bp8==sb8(8) )
//      ppf=pbf(1)
//      pfa(l+1)=int((ppf-1)/23)+1; bfa(l+1)=ppf-pfa(l+1)*23+23
//// note that this also works since 48 and 23 happen to be prime
////  among themselves
//end
//disp(int8(pfa))
//disp(int8(bfa))
//
/// The explanation is: for all the 256 possibilities of having
/// each of the 8 vertices of the cell above/below the threshold,
///  specify which is the relevant triangle pattern (pf) and 
/// arrangement (bf)

pfa=[1 1 7 1 3 1 3 1 19 13 7 7 15 13 9 1 20 14 8 8 3 ..
  20 3 1 25 25 9 7 15 13 9 36 4 2 4 2 5 13 5 1 4 ..
  19 4 2 21 13 5 40 16 14 10 2 22 14 6 39 16 14 10 ..
  35 5 47 41 41 23 1 11 1 17 23 11 3 27 13 7 7 27 15 ..
  9 34 28 14 8 8 17 1 11 37 25 25 7 5 27 10 2 34 18 ..
  24 12 4 17 17 5 37 18 2 12 38 17 44 32 38 28 16 10 ..
  33 18 43 31 37 28 9 1 33 7 7 31 43 43 31 7 7 33 1 ..
  9 1 37 31 43 19 33 13 21 28 38 32 44 20 38 32 3 18 ..
  37 25 37 17 37 12 24 18 34 2 10 2 5 1 5 25 37 31 4 ..
  17 33 8 14 28 34 14 22 27 34 7 13 27 38 11 23 17 1 ..
  11 1 23 41 41 47 1 35 35 23 16 39 39 39 22 27 10 14 ..
  16 40 40 40 21 42 4 19 4 37 5 13 5 2 4 2 4 36 36 ..
  24 15 17 9 25 25 41 3 20 3 8 8 14 20 28 9 13 15 7 ..
  7 13 19 1 3 1 3 1 7 1 1 ] 
 
bfa=[1 2 2 3 2 4 3 5 2 3 4 5 3 5 5 9 2 3 4 5 6 7 7 ..
  13 4 5 8 14 7 13 15 19 2 4 3 5 4 8 5 14 6 7 7 13 ..
  7 15 13 19 3 5 5 9 7 15 13 19 7 13 15 19 11 17 17 ..
  21 2 6 4 7 3 7 5 13 4 7 8 15 5 13 14 19 4 7 8 15 ..
  7 11 15 17 8 15 12 16 15 17 16 20 3 7 5 13 5 15 9 ..
  19 7 11 15 17 13 17 19 21 5 13 14 19 13 17 19 21 15 ..
  17 16 20 17 18 20 22 2 4 6 7 4 8 7 15 3 5 7 13 5 ..
  14 13 19 3 5 7 13 7 15 10 17 5 9 15 19 13 19 17 21 ..
  4 8 7 15 8 12 15 16 7 15 10 17 15 16 17 20 5 14 13 ..
  19 15 16 17 20 13 19 17 21 17 20 18 22 3 7 7 10 5 ..
  15 13 17 5 13 15 17 9 19 19 21 5 13 15 17 13 17 17 ..
  18 14 19 16 20 19 21 20 22 5 15 13 17 14 16 19 20 ..
  13 17 17 18 19 20 21 22 9 19 19 21 19 20 21 22 19 ..
  21 20 22 21 22 22 23 ] 

// also this sits better outside of the loops:
v=ones(255,24);
for l=(2:255)
  v(l,:)=ppinv(pfa(l),vtrian(bfa(l)*ppar(pfa(l))+23,:));
end

nfta=nft(bfa(1:256));

// now scan every cube of the 3d hypermat

nf=0

// a few definitions brought out of the loops
u1=int((0:35)/3 +1)*2-1; u2=u1+1;
v1=1:2:23; 
//v2=v1+1; v12=[v1,v2]; va=1:12; vb=va+12;
uu1=4*u1; uu2=4*u2; 
uv2=[v1*4+1;v1*4+2;v1*4+3]; uv1=uv2-4;

xyz=[];
// to dimension xyz maximally (e.g. (9,4*(nxm)*(nym)*(nzm)) ) 
// from start would be a waste of memory, but to grow it for
// every added triangle would be quite slow. I think that a good 
// compromise is to resize xx,yy,zz at each step of the i loop
// (there xyz needs to grow only from (9,nt-1) to (9,4*nxm*nym+nt-1) )


nt=0

if datastruct=="h" then
  nxm=nx-1; nym=ny-1; nzm=nz-1
  nx1=1:nxm; ny1=1:nym; nz1=1:nzm; ny2=ny1+1; nz2=nz1+1;
  u=1:3; nyzm=nym*nzm; Ny=1:ny; Nz=1:nz; js=(0:7)*nym;
  yy8=[y(ny1);y(ny1);y(ny2);y(ny2);y(ny1);y(ny1);y(ny2);y(ny2)]'
  zz8=[z(nz1);z(nz1);z(nz1);z(nz1);z(nz2);z(nz2);z(nz2);z(nz2)]'
// sort the points of the cube in canonical order - the
//  outer in the loop the faster (?)
  for i=nx1
// I scan the hypermatrix with i outer, try to collapse j and k
   i1=i+1;
   x8=x([i,i1,i1,i,i,i1,i1,i])
//   x8=x([0 1 1 0 0 1 1 0]+i)
   s8i=matrix(s(i,Ny,Nz),ny,nz); //s8i=s8i(:,:)
   s8i1=matrix(s(i1,Ny,Nz),ny,nz); //s8i1=s8i1(:,:)
// s8i1(:,:) etc instead of s8i1 for compatibility with scilab 2.4.1
   ss8=double([s8i(ny1,nz1); s8i1(ny1,nz1); s8i1(ny2,nz1); s8i(ny2,nz1); ..
             s8i(ny1,nz2); s8i1(ny1,nz2); s8i1(ny2,nz2); s8i(ny2,nz2)])'-s0; 
// this matrix ss8 now contains all the necessary values of s in the
// planes i and i+1. I play index tricks, and try to minimize the
// use of hypermatrices, as calculations with them are considerably
// slower. 
   ll=matrix(matrix(bool2s(ss8>0),nyzm,8)*pow2+1,nzm,nym);
// Recognize which pattern of >0 and <0 appears in each cell of 
// this plane: boolean((ny-1)*(nz-1),8) --> integers in a snap, 
// for the whole plane.
// Usually most of the cells have all 8 values either >s0 or <s0,
// and don't need to be considered for isosurfacelets. So the following
//  loop is only on those cells which have some vertices on 
//  opposite sides of the s0 divider.
   [kk,jj]=find(ll>1 & ll<256); njk=length(jj); 
// making space in xyz
   xyz=[xyz(1:9*nt);zeros(36*njk,1)];
   for q=1:njk
       j=jj(q); k=kk(q);
       l=ll(k,j)
       nftbf=nfta(l);
// pick up the relevant vertices for constructing triangles
// and generate the triangular facelets
       X=[x8;yy8(j,:);zz8(k,:);ss8(k,js+j)]
//   X 4x8 stores {x,y,z,s-s0} of the 8 vertices of the cell.
//timer()
       Y=X(:,v(l,:))  // Y(4,24)
//the following criptic lines interpolate linearly along the relevant
// sides, and return nftbf triples of x,y,z, vertices of the triangles.
// Matrix flattening is used to increase speed. Index arrays like
// u,u1,u2,v1,v2,uu1,uu2,uv1,uv2, have been defined above.
       S1=Y(uu1); S2=Y(uu2);
       D=S2-S1; N=S2.*Y(uv1)-S1.*Y(uv2);
// denominators are never zero if Y is formed with the right points
//  and if only the first nftbf triangles are computed
       uu=1:9*nftbf;
       xyz(9*nt+uu)=N(uu)./D(uu)
//t=t+timer()
// There is a redundancy of operations, as points are
// always aligned (have two equal coordinates), so will be their
// mean point, and there would be no need to compute it. I'm afraid 
// that any catch of that will only overburden the snippet.
       nt=nt+nftbf
   end
  end
else
  ll=matrix(bool2s(field(kvert)>f0),ncell,8)*pow2+1
  jj=find(ll>1 & ll<256);
  xyz=zeros(9*sum(nfta(ll(jj))),1);
  for j=jj
       l=ll(j); k=kvert(j,:); nftbf=nfta(l);
       X=[dcorvg(k,1)';dcorvg(k,2)';dcorvg(k,3)';field(k)'-f0]
//   X 4x8 stores {x,y,z,s-s0} of the 8 vertices of the cell.
       Y=X(:,v(l,:))  // Y(4,24)
       S1=Y(uu1); S2=Y(uu2);
       D=S2-S1; N=S2.*Y(uv1)-S1.*Y(uv2);
       uu=1:9*nftbf;
       xyz(9*nt+uu)=N(uu)./D(uu)
       nt=nt+nftbf
  end
end
// output only the xx,yy,zz needed
u=1:nt
xyz=matrix(xyz,9,-1)
xx=xyz(1:3:7,u); yy=xyz(2:3:8,u); zz=xyz(3:3:9,u); 


//write(%io(2),'time spent in the selected routine:'+string(t))

endfunction

function shadesurf(x,y,z,a,b,theta,alpha,leg,flag,ebox)

//  plot of a surface with lighting due to an infinitely
//    distant light source (the sun), no shading
//
// usage: shadesurf(x,y,z,[a,b,theta,alpha,leg,flag,ebox])
//
//     all arguments like those of plot3d1 (in fact, that's
//     a little more than a wrapper), except for:
//     a,b:  director cosines of the sun direction
//   


[lhs,rhs]=argn(0);

if rhs==0 then
//demo
  disp "Demo of shadesurf:"
  disp "  x=linspace(-2.5,2.5,50);"
  disp "  a=exp(-(x''*x).^2);"
  disp "  xbasc(); shadesurf(x,x,a,-1,1,65,60);"
  x=linspace(-2.5,2.5,35); 
  a=exp(-(x'*x).^2);
  xbasc(); shadesurf(x,x,a,-1,1,65,60); 
  return
end
if rhs<4 then a=0; end
if rhs<5 then b=.2; end
if rhs<6 then theta=45; end
if rhs<7 then alpha=35; end
if rhs<8 then leg='X@Y@Z'; end
if rhs<9 then flag=[-1 2 4]'; end
if rhs<10 then ebox=[0 0 0 1 1 1]'; end


nc=xget("lastpattern");
// beware - drivers other than 'REC' have their bugs about colormaps


if size(x,1)==1 | size(x,2)==1 then
// style like plot3d(x,y,z): x(1:nx), y(1:ny), z(1:nx,1:ny)
   [xx,yy,zz]=genfac3d(x,y,z);
   nn1=(zz(1,:)+zz(2,:)-zz(3,:)-zz(4,:))./(xx(4,:)-xx(1,:));
   nn2=(zz(1,:)-zz(2,:)-zz(3,:)+zz(4,:))./(yy(2,:)-yy(1,:));
   nn3=1
else
// style like plot3d(xx,yy,zz): facelets (thought primarily for
// triangles)
  xx=x; yy=y; zz=z;
  nn1=-((yy(1,:)-yy(2,:)).*(zz(2,:)-zz(3,:))-...
      (yy(2,:)-yy(3,:)).*(zz(1,:)-zz(2,:)));
  nn2=((xx(1,:)-xx(2,:)).*(zz(2,:)-zz(3,:))-...
      (xx(2,:)-xx(3,:)).*(zz(1,:)-zz(2,:)));
  nn3=-((xx(1,:)-xx(2,:)).*(yy(2,:)-yy(3,:))-...
      (xx(2,:)-xx(3,:)).*(yy(1,:)-yy(2,:)));
  if size(xx,1)==3 then
// so far scilab (up to 2.5) has the bug of the extra line of the
// 3d triangle
    xx=[xx;xx(1,:)]; yy=[yy;yy(1,:)]; zz=[zz;zz(1,:)]; 
  end
end

degfac=find(abs(nn1)+abs(nn2)+abs(nn3)<10*%eps);
// this can happen for instance if there are coincident vertices;
// it shouldn't, if xx, yy, zz are constructed properly, but
// I want shadesurf to get through it anyway
nn3(degfac)=1

//cc=sqrt(nn1.^2+nn2.^2+nn3.^2);
// color proportional to the slope

//cc=1../max(a*nn1+b*nn2,1);
// color inversely proportional to the lightened side

cc=(a*nn1+b*nn2+nn3)./sqrt(nn1.^2+nn2.^2+nn3.^2) ./sqrt(a^2+b^2+1);
// color proportional to the sine of the angle (local normal)^(sun)
// the direction of the sun is [a, b, 1]

if max(cc)~=min(cc) then
   cc=(-min(cc)+cc)/(max(cc)-min(cc))*(nc-1)+1;
else
// this can happen if the surface is a plane, and I want this to
//  be handled
   cc=ones(size(cc,1),size(cc,2))*nc;
end

plot3d1(xx,yy,list(zz,cc),theta,alpha,leg,flag,ebox)

// to draw also the external contour of the surface would be
//  nice, but requires more tracking of the external edges

endfunction