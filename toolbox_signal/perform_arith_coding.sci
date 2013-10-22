function [y,nbr_bits] = perform_arith_coding(xC, dir)

// perform_arith_coding - perform adaptive arithmetic coding
//
// 	[y,nbr_bits] = perform_arith_coding(x, dir);
//
//   dir=1 for encoding, dir=-1 for decoding.
//
// Interface to the code of (c) Karl Skretting
//
//   Copyright (c) 2008 Gabriel Peyre

//
// ------------------------------------------------------------------
// Arguments:
//  y        a column vector of non-negative integers (bytes) representing
//           the code, 0 <= y(i) <= 255.
//  Res      a matrix that sum up the results, size is (NumOfX+1)x4
//           one line for each of the input sequences, the columns are
//           Res(:,1) - number of elements in the sequence
//           Res(:,2) - unused (=0)
//           Res(:,3) - bits needed to code the sequence
//           Res(:,4) - bit rate for the sequence, Res(:,3)/Res(:,1)
//           Then the last line is total (which include bits needed to store NumOfX)
//  x       a cell array of column vectors of integers representing the
//           symbol sequences. (should not be to large integers)
//           If only one sequence is to be coded, we must make the cell array
//           like: xC=cell(2,1); xC{1}=x; // where x is the sequence
// ------------------------------------------------------------------
// Note: this routine is extremely slow since it is all Matlab code

// SOME NOTES ON THE FUNCTION
// This function is almost like Arith06, but some important changes have
// been done. Arith06 is buildt almost like Huff06, but this close connection
// is removed in Arith07. This imply that to understand the way Arith06
// works you should read the dokumentation for Huff06 and especially the
// article on Recursive Huffman Coding. To understand how Arith07 works it is
// only confusing to read about the recursive Huffman coder, Huff06.
//
// Arith07 code each of the input sequences in xC in sequence, the
// method used for each sequence depends on what kind (xType) the
// sequence is. We have
//  xType  Explanation
//   0     Empty sequence (L=0)
//   1     Only one symbol (L=1) and x>1
//   9     Only one symbol (L=1) and x=0/1
//  11     Only one symbol (L=1) and x<0
//   2     all symbols are equal, L>1, x(1)=x(2)=...=x(L)>1
//   8     all symbols are equal, L>1, x(1)=x(2)=...=x(L)=0/1
//  10     all symbols are equal, L>1, x(1)=x(2)=...=x(L)<0
//   3     non-negative integers, 0<=x(l)<=1000
//   4     not to large integers, -1000<=x(l)<=1000
//   5     large non-negative integers possible, 0<=x(l)
//   6     also possible with large negative numbers
//   7     A binary sequence, x(l)=0/1.
// Many functions are defined in this m-file:
//  PutBit, GetBit     read and write bit from/to bitstream
//  PutABit, GetABit   Arithmetic encoding of a bit, P{0}=P{1}=0.5
//  PutVLIC, GetVLIC   Variable Length Integer Code
//  PutS(x,S,C), GetS(S,C)  A symbol x, which is in S, is aritmetic coded
//                     according to the counts given by C.
//                     Ex: A binary variable with the givene probabilities
//                     P{0}=0.8 and P{1}=0.2, PutS(x,[0,1],[4,1]).
//  PutN(x,N), GetN(N) Arithmetic coding of a number x in range 0<=x<=N where
//                     we assume equal probability, P{0}=P{1}=...=P{N}=1/(N+1)


xC = uint32(xC);

if dir==1
//    xC = {xC};
end

//----------------------------------------------------------------------
// Copyright (c) 1999-2001.  Karl Skretting.  All rights reserved.
// Hogskolen in Stavanger (Stavanger University), Signal Processing Group
// Mail:  karl.skretting@tn.his.no   Homepage:  http://www.ux.his.no/~karlsk/
//
// HISTORY:
// Ver. 1.0  19.05.2001  KS: Function made (based on Arith06 and Huff06)
//----------------------------------------------------------------------

// these global variables are used to read from or write to the compressed sequence
global y Byte BitPos
// and these are used by the subfunctions for arithmetic coding
global high low range ub hc lc sc K code

Mfile='Arith07';
K=24;  // number of bits to use in integers  (4 <= K <= 24)
Display=0;      // display progress and/or results

// check input and output arguments, and assign values to arguments
if (argn(2) < 1);
    error('function must have input arguments, see help.');
end
if (argn(1) < 1);
    error('function must have output arguments, see help.');
end

if dir==0 // (~iscell(xC))
    Encode=0;Decode=1;
    y=xC(:);            // first argument is y
    y=[y;0;0;0;0];      // add some zeros to always have bits available
else
    Encode=1;Decode=0;
    NumOfX = 1; // length(xC);
end

Byte=0;BitPos=1;         // ready to read/write into first position
low=0;high=2^K-1;ub=0;   // initialize the coder
code=0;
// we just select some probabilities which we belive will work quite well
TypeS=[  3,  4,  5,  7,  6,  1,  9, 11,  2,  8, 10,  0];
TypeC=[100, 50, 50, 50, 20, 10, 10, 10,  4,  4,  2,  1];

if Encode
    Res=zeros(NumOfX,4);
    // initalize the global variables
    y=zeros(10,1);    // put some zeros into y initially
    // start encoding, first write VLIC to give number of sequences
    PutVLIC(NumOfX);
    // now encode each sequence continuously
    Ltot=0;
    for num=1:NumOfX
        x=xC(:); // {num};
//        x=full(x(:));        // make sure x is a non-sparse column vector
        L=length(x);
        Ltot=Ltot+L;
        y=[y(1:Byte);zeros(50+2*L,1)];  // make more space available in y
        StartPos=Byte*8-BitPos+ub;      // used for counting bits
        // find what kind of sequenqe we have, save this in xType
        if (L==0)
            xType=0;
        elseif L==1
            if (x==0) | (x==1)       // and it is 0/1
                xType=9;
            elseif (x>1)             // and it is not 0/1 (and positive)
                xType=1;
            else                     // and it is not 0/1 (and negative)
                xType=11;
            end
        else
            // now find some more info about x
            maxx=max(x);
            minx=min(x);
            rangex=maxx-minx+1;
            if (rangex==1)                    // only one symbol
                if (maxx==0) | (maxx==1)       // and it is 0/1
                    xType=8;
                elseif (maxx>1)                // and it is not 0/1 (and positive)
                    xType=2;
                else                           // and it is not 0/1 (and negative)
                    xType=10;
                end
            elseif (minx == 0) & (maxx == 1)   // a binary sequence
                xType=7;
            elseif (minx >= 0) & (maxx <= 1000)
                xType=3;
            elseif (minx >= 0)
                xType=5;
            elseif (minx >= -1000) & (maxx <= 1000)
                xType=4;
            else
                xType=6;
            end
        end        // if (L==0)
        //
        if sum(xType==[4,6])        // negative values are present
            I=find(x);               // non-zero entries in x
            Sg=(sign(x(I))+1)/2;     // the signs will be needed later, 0/1
            x=abs(x);
        end
        if sum(xType==[5,6])        // we take the logarithms of the values
            I=find(x);               // non-zero entries in x
            xa=x;                    // additional bits
            x(I)=floor(log2(x(I)));
            xa(I)=xa(I)-2.^x(I);
            x(I)=x(I)+1;
        end
        //  now do the coding of this sequence
        PutS(xType,TypeS,TypeC);
        if (xType==1)       // one symbol and x=x(1)>1
            PutVLIC(x-2);
        elseif (xType==2)   // L>1 but only one symbol, x(1)=x(2)=...=x(L)>1
            PutVLIC(L-2);
            PutVLIC(x(1)-2);
        elseif sum(xType==[3,4,5,6])  // now 'normalized' sequences: 0 <= x(i) <= 1000
            PutVLIC(L-2);
            M=max(x);
            PutN(M,1000);         // some bits for M
            // initialize model
            T=[ones(M+2,1);0];
            Tu=flipud((-1:(M+1))');   // (-1) since ESC never is used in Tu context
            // and code the symbols in the sequence x
            for l=1:L
                sc=T(1);
                m=x(l);
                hc=T(m+1);lc=T(m+2);
                if hc==lc      // unused symbol, code ESC symbol first
                    hc=T(M+2);lc=T(M+3);
                    EncodeSymbol;      // code escape with T table
                    sc=Tu(1);hc=Tu(m+1);lc=Tu(m+2);  // symbol with Tu table
                    Tu(1:(m+1))=Tu(1:(m+1))-1;       // update Tu table
                end
                EncodeSymbol;  // code actual symbol with T table (or Tu table)
                // update T table, MUST be identical in Encode and Decode
                // this avoid very large values  in T even if sequence is very long
                T(1:(m+1))=T(1:(m+1))+1;
                if (rem(l,5000)==0)
                    dT=T(1:(M+2))-T(2:(M+3));
                    dT=floor(dT*7/8+1/8);
                    for m=(M+2):(-1):1; T(m)=T(m+1)+dT(m); end;
                end
            end          // for l=1:L
            // this end the "elseif sum(xType==[3,4,5,6])"-clause
        elseif (xType==7)   // L>1 and   0 <= x(i) <= 1
            PutVLIC(L-2);
            EncodeBin(x,L);    // code this sequence a special way
        elseif (xType==8)   // L>1 and   0 <= x(1)=x(2)=...=x(L) <= 1
            PutVLIC(L-2);
            PutABit(x(1));
        elseif (xType==9)   // L=1 and   0 <= x(1) <= 1
            PutABit(x(1));
        elseif (xType==10)       // L>1 and  x(1)=x(2)=...=x(L) <= -1
            PutVLIC(L-2);
            PutVLIC(-1-x(1));
        elseif (xType==11)       // L=1 and   x(1) <= -1
            PutVLIC(-1-x);
        end         // if (xType==1)
        // additional information should be coded as well
        if 0         // first the way it is not done any more
            if sum(xType==[4,6])           // sign must be stored
                for i=1:length(Sg); PutABit(Sg(i)); end;
            end
            if sum(xType==[5,6])        // additional bits must be stored
                for i=1:L
                    for ii=(x(i)-1):(-1):1
                        PutABit(bitget(xa(i),ii));
                    end
                end
            end
        else         // this is how we do it
            if sum(xType==[4,6])           // sign must be stored
                EncodeBin(Sg,length(I));    // since length(I)=length(Sg)
            end
            if sum(xType==[5,6])        // additional bits must be stored
                b=zeros(sum(x)-length(I),1);  // number of additional bits
                bi=0;
                for i=1:L
                    for ii=(x(i)-1):(-1):1
                        bi=bi+1;
                        b(bi)=bitget(xa(i),ii);
                    end
                end
                if (bi~=(sum(x)-length(I)))
                    error([Mfile,': logical error, bi~=(sum(x)-length(I)).']);
                end
                EncodeBin(b,bi);    // since bi=(sum(x)-length(I))
            end
        end
        //
        EndPos=Byte*8-BitPos+ub;    // used for counting bits
        bits=EndPos-StartPos;
        Res(num,1)=L;
        Res(num,2)=0;
        Res(num,3)=bits;
        if L>0; Res(num,4)=bits/L; else Res(num,4)=bits; end;
    end         // for num=1:NumOfX
    // flush the arithmetic coder
    PutBit(bitget(low,K-1));
    ub=ub+1;
    while ub>0
        PutBit(~bitget(low,K-1));
        ub=ub-1;
    end
    // flush is finished
    y=y(1:Byte);
    // varargout(1) = {y};
    if (nargout >= 2)
        // now calculate results for the total
        Res(NumOfX+1,1)=Ltot;
        Res(NumOfX+1,2)=0;
        Res(NumOfX+1,3)=Byte*8;
        if (Ltot>0); Res(NumOfX+1,4)=Byte*8/Ltot; else Res(NumOfX+1,4)=Byte*8; end;
        // varargout(2) = {Res};
    end
end         // if Encode

if Decode
    for k=1:K
        code=code*2;
        code=code+GetBit;   // read bits into code
    end
    NumOfX=GetVLIC;   // first read number of sequences
    xC=cell(NumOfX,1);
    for num=1:NumOfX
        // find what kind of sequenqe we have, xType, stored first in sequence
        xType=GetS(TypeS,TypeC);
        // now decode the different kind of sequences, each the way it was stored
        if (xType==0)       // empty sequence, no more symbols coded
            x=[];
        elseif (xType==1)       // one symbol and x=x(1)>1
            x=GetVLIC+2;
        elseif (xType==2)   // L>1 but only one symbol, x(1)=x(2)=...=x(L)>1
            L=GetVLIC+2;
            x=ones(L,1)*(GetVLIC+2);
        elseif sum(xType==[3,4,5,6])  // now 'normalized' sequences: 0 <= x(i) <= 1000
            L=GetVLIC+2;
            x=zeros(L,1);
            M=GetN(1000);      // M is max(x)
            // initialize model
            T=[ones(M+2,1);0];
            Tu=flipud((-1:(M+1))');   // (-1) since ESC never is used in Tu context
            // and decode the symbols in the sequence x
            for l=1:L
                sc=T(1);
                range=high-low+1;
                counts=floor(( (code-low+1)*sc-1 )/range);
                m=2; while (T(m)>counts); m=m+1; end;
                hc=T(m-1);lc=T(m);m=m-2;
                RemoveSymbol;
                if (m>M)     // decoded ESC symbol, find symbol from Tu table
                    sc=Tu(1);range=high-low+1;
                    counts=floor(( (code-low+1)*sc-1 )/range);
                    m=2; while (Tu(m)>counts); m=m+1; end;
                    hc=Tu(m-1);lc=Tu(m);m=m-2;
                    RemoveSymbol;
                    Tu(1:(m+1))=Tu(1:(m+1))-1;   // update Tu table
                end
                x(l)=m;
                // update T table, MUST be identical in Encode and Decode
                // this avoid very large values  in T even if sequence is very long
                T(1:(m+1))=T(1:(m+1))+1;
                if (rem(l,5000)==0)
                    dT=T(1:(M+2))-T(2:(M+3));
                    dT=floor(dT*7/8+1/8);
                    for m=(M+2):(-1):1; T(m)=T(m+1)+dT(m); end;
                end
            end          // for l=1:L
            // this end the "elseif sum(xType==[3,4,5,6])"-clause
        elseif (xType==7)   // L>1 and   0 <= x(i) <= 1
            L=GetVLIC+2;
            x=DecodeBin(L);    // decode this sequence a special way
        elseif (xType==8)   // L>1 and   0 <= x(1)=x(2)=...=x(L) <= 1
            L=GetVLIC+2;
            x=ones(L,1)*GetABit;
        elseif (xType==9)   // L=1 and   0 <= x(1) <= 1
            x=GetABit;
        elseif (xType==10)       // L>1 and  x(1)=x(2)=...=x(L) <= -1
            L=GetVLIC+2;
            x=ones(L,1)*(-1-GetVLIC);
        elseif (xType==11)       // L=1 and   x(1) <= -1
            x=(-1-GetVLIC);
        end         // if (xType==0)
        // additional information should be decoded as well
        L=length(x);
        I=find(x);
        if 0
            if sum(xType==[4,6])           // sign must be retrieved
                Sg=zeros(size(I));
                for i=1:length(I); Sg(i)=GetABit; end;   // and the signs   (0/1)
                Sg=Sg*2-1;                               // (-1/1)
            end
            if sum(xType==[5,6])        // additional bits must be retrieved
                xa=zeros(L,1);
                for i=1:L
                    for ii=2:x(i)
                        xa(i)=2*xa(i)+GetABit;
                    end
                end
                x(I)=2.^(x(I)-1);
                x=x+xa;
            end
        else
            if sum(xType==[4,6])           // sign must be retrieved
                Sg=DecodeBin(length(I));    // since length(I)=length(Sg)
                Sg=Sg*2-1;                  // (-1/1)
            end
            if sum(xType==[5,6])        // additional bits must be retrieved
                bi=sum(x)-length(I);     // number of additional bits
                b=DecodeBin(bi);
                bi=0;
                xa=zeros(L,1);
                for i=1:L
                    for ii=2:x(i)
                        bi=bi+1;
                        xa(i)=2*xa(i)+b(bi);
                    end
                end
                x(I)=2.^(x(I)-1);
                x=x+xa;
            end
        end
        if sum(xType==[4,6])           // sign must be used
            x(I)=x(I).*Sg;
        end
        // now x is the retrieved sequence
        xC = x; // {num}=x;
    end         // for num=1:NumOfX
    // varargout(1) = {xC};
end

if dir==1
    nbr_bits = Res(1,3);
else
    // undef for decoding
    nbr_bits = -1;
    y = xC; // {1};
end


endfunction


//----------------------------------------------------------------------
//----------------------------------------------------------------------

// --- The functions for binary sequences: EncodeBin and DecodeBin -------
// These function may call themselves recursively
function EncodeBin(x,L)
global y Byte BitPos
global high low range ub hc lc sc K code

Display=0;
x=x(:);
if (length(x)~=L); error('EncodeBin: length(x) not equal L.'); end;

// first we try some different coding methods to find out which one
// that might do best. Many more methods could have been tried, for
// example methods to check if x is a 'byte' sequence, or if the bits
// are grouped in other ways. The calling application is the best place
// to detect such dependencies, since it will (might) know the process and
// possible also its statistical (or deterministic) properties.
// Here we just check some few coding methods: direct, split, diff, diff+split
// The main variables used for the different methods are:
//  direct: x, I, J, L11, b0
//  split: x is split into x1 and x2, L1, L2, b1 is bits needed
//  diff: x3 is generated from x, I3, J3 L31, b2
//  diff+split: x3 is split into x4 and x5, L4, L5, b3
//
MetS=[0,1,2,3];       // the different methods, direct, split, diff, diff+split
MetC=[9,3,3,1];       // and the counts (which gives the probabilities)
// first set how many bits needed to code the method
b0=log2(sum(MetC))-log2(MetC(1));
b1=log2(sum(MetC))-log2(MetC(2));
b2=log2(sum(MetC))-log2(MetC(3));
b3=log2(sum(MetC))-log2(MetC(4));

I=find(x(1:(L-1))==1);      // except last element in x
J=find(x(1:(L-1))==0);
L11=length(I);
// perhaps is 'entropy' interesting
// N1=L11+x(L);
// N0=L-N1;
// b=N1*log2(L/N1)+N0*log2(L/N0);
// disp(['L=',int2str(L),',  N1=',int2str(N1),',  (e)bits=',num2str(b)]);
b0=b0+BitEst(L,L11+x(L));   // bits needed for the sequence without splitting
if Display
    disp(['EncodeBin: x of length ',int2str(L),' and ',int2str(L11+x(L)),...
        ' ones (p=',num2str((L11+x(L))/L,'//1.3f'),...
        ') can be coded by ',num2str(b0,'//6.0f'),' bits.']);
end
// diff with a binary sequence is to indicate wether or not a symbol is
// the same as preceeding symbol or not, x(0) is assumed to be '0' (zero).
// This is the DPCM coding scheme on a binary sequence
x3=abs(x-[0;x(1:(L-1))]);
I3=find(x3(1:(L-1))==1);      // except last element in x3
J3=find(x3(1:(L-1))==0);
L31=length(I3);
b2=b2+BitEst(L,L31+x3(L));  // bits needed for the sequence without splitting
if Display
    disp(['EncodeBin: diff x, L=',int2str(L),' gives ',int2str(L31+x3(L)),...
        ' ones (p=',num2str((L31+x3(L))/L,'//1.3f'),...
        ') can be coded by ',num2str(b2,'//6.0f'),' bits.']);
end
//
if (L>40)
    // only now we try to split the sequences x and x3
    if (L11>3) & ((L-L11)>3)
        // try to split x into x1 and x2, depending on previous symbol
        if L11<(L/2)
            x1=x(I+1);
            x2=[x(1);x(J+1)];
            L1=L11;L2=L-L11;
        else
            x1=[x(1);x(I+1)];
            x2=x(J+1);
            L1=L11+1;L2=L-L11-1;
        end
        b11=BitEst(L1,length(find(x1)));  // bits needed for x1
        b12=BitEst(L2,length(find(x2)));  // bits needed for x2
        // b1 is bits to code: Method, L11, x1 and x2
        b1=b1+log2(L)+b11+b12;
        if Display
            disp(['EncodeBin, x -> x1+x2:      lengths are ',int2str(L1),'+',int2str(L2),...
                '  bits are ',num2str(b11,'//6.0f'),'+',num2str(b12,'//6.0f'),...
                '.  Total is ',num2str(b1,'//6.0f'),' bits.']);
        end
    else
        b1=b0+1;      // just to make this larger
    end
    if (L31>3) & ((L-L31)>3)
        // try to split x3 into x4 and x5, depending on previous symbol
        if L31<(L/2)
            x4=x3(I3+1);
            x5=[x3(1);x3(J3+1)];
            L4=L31;L5=L-L31;
        else
            x4=[x3(1);x3(I3+1)];
            x5=x3(J3+1);
            L4=L31+1;L5=L-L31-1;
        end
        b31=BitEst(L4,length(find(x4)));  // bits needed for x4
        b32=BitEst(L5,length(find(x5)));  // bits needed for x5
        // b3 is bits to code: Method, L31, x4 and x5
        b3=b3+log2(L)+b31+b32;
        if Display
            disp(['EncodeBin, diff x -> x4+x5: lengths are ',int2str(L4),'+',int2str(L5),...
                '  bits are ',num2str(b31,'//6.0f'),'+',num2str(b32,'//6.0f'),...
                '.  Total is ',num2str(b3,'//6.0f'),' bits.']);
        end
    else
        b3=b2+1;      // just to make this larger
    end
else
    b1=b0+1;      // just to make this larger
    b3=b2+1;      // just to make this larger
end
// now code x by the best method of those investigated
[b,MetI]=min([b0,b1,b2,b3]);
MetI=MetI-1;
PutS(MetI,MetS,MetC);   // code which method to use
//
if MetI==0
    // code the sequence x
    N1=L11+x(L);
    N0=L-N1;
    PutN(N1,L);            // code N1, 0<=N1<=L
    for n=1:L
        if ~(N0*N1); break; end;
        PutS(x(n),[0,1],[N0,N1]);   // code x(n)
        N0=N0-1+x(n);N1=N1-x(n);    // update model (of rest of the sequence)
    end
elseif MetI==1
    // code x1 and x2
    clear x4 x5 x3 x I3 J3 I J
    PutN(L11,L-1);         // 0<=L11<=(L-1)
    EncodeBin(x1,L1);
    EncodeBin(x2,L2);
elseif MetI==2
    // code the sequence x3
    N1=L31+x3(L);
    N0=L-N1;
    PutN(N1,L);            // code N1, 0<=N1<=L
    for n=1:L
        if ~(N0*N1); break; end;
        PutS(x3(n),[0,1],[N0,N1]);   // code x3(n)
        N0=N0-1+x3(n);N1=N1-x3(n);    // update model (of rest of the sequence)
    end
elseif MetI==3
    // code x4 and x5
    clear x1 x2 x3 x I3 J3 I J
    PutN(L31,L-1);         // 0<=L31<=(L-1)
    EncodeBin(x4,L4);
    EncodeBin(x5,L5);
end
return    // end of EncodeBin

endfunction

function x = DecodeBin(L)
global y Byte BitPos
global high low range ub hc lc sc K code

// these must be as in EncodeBin
MetS=[0,1,2,3];       // the different methods, direct, split, diff, diff+split
MetC=[9,3,3,1];       // and the counts (which gives the probabilities)
MetI=GetS(MetS,MetC);    // encode which method to use

if (MetI==1) | (MetI==3)   // a split was done
    L11=GetN(L-1);         // 0<=L11<=(L-1)
    if L11<(L/2)
        L1=L11;L2=L-L11;
    else
        L1=L11+1;L2=L-L11-1;
    end
    x1=DecodeBin(L1);
    x2=DecodeBin(L2);
    // build sequence x from x1 and x2
    x=zeros(L,1);
    if L11<(L/2)
        x(1)=x2(1);
        n1=0;n2=1;   // index for the last in x1 and x2
    else
        x(1)=x1(1);
        n1=1;n2=0;   // index for the last in x1 and x2
    end
    for n=2:L
        if (x(n-1))
            n1=n1+1;
            x(n)=x1(n1);
        else
            n2=n2+1;
            x(n)=x2(n2);
        end
    end
else                      // no split
    N1=GetN(L);
    N0=L-N1;
    x=zeros(L,1);
    for n=1:L
        if (N0==0); x(n:L)=1; break; end;
        if (N1==0); break; end;
        x(n)=GetS([0,1],[N0,N1]);   // decode x(n)
        N0=N0-1+x(n);N1=N1-x(n);    // update model (of rest of the sequence)
    end
end

if (MetI==2) | (MetI==3)   // x is diff coded
    for n=2:L
        x(n)=x(n-1)+x(n);
    end
    x=rem(x,2);
end

return    // end of DecodeBin

endfunction

// ------- Other subroutines ------------------------------------------------

// Functions to write and read a Variable Length Integer Code word
// This is a way of coding non-negative integers that uses fewer
// bits for small integers than for large ones. The scheme is:
//   '00'   +  4 bit  - integers from 0 to 15
//   '01'   +  8 bit  - integers from 16 to 271
//   '10'   + 12 bit  - integers from 272 to 4367
//   '110'  + 16 bit  - integers from 4368 to 69903
//   '1110' + 20 bit  - integers from 69940 to 1118479
//   '1111' + 24 bit  - integers from 1118480 to 17895695
//   not supported  - integers >= 17895696 (=2^4+2^8+2^12+2^16+2^20+2^24)
function PutVLIC(N)
global y Byte BitPos
global high low range ub hc lc sc K code
if (N<0)
    error('Arith07-PutVLIC: Number is negative.');
elseif (N<16)
    PutABit(0);PutABit(0);
    for (i=4:-1:1); PutABit(bitget(N,i)); end;
elseif (N<272)
    PutABit(0);PutABit(1);
    N=N-16;
    for (i=8:-1:1); PutABit(bitget(N,i)); end;
elseif (N<4368)
    PutABit(1);PutABit(0);
    N=N-272;
    for (i=12:-1:1); PutABit(bitget(N,i)); end;
elseif (N<69940)
    PutABit(1);PutABit(1);PutABit(0);
    N=N-4368;
    for (i=16:-1:1); PutABit(bitget(N,i)); end;
elseif (N<1118480)
    PutABit(1);PutABit(1);PutABit(1);PutABit(0);
    N=N-69940;
    for (i=20:-1:1); PutABit(bitget(N,i)); end;
elseif (N<17895696)
    PutABit(1);PutABit(1);PutABit(1);PutABit(1);
    N=N-1118480;
    for (i=24:-1:1); PutABit(bitget(N,i)); end;
else
    error('Arith07-PutVLIC: Number is too large.');
end

endfunction




// Aritmetic coding of a number or symbol x, the different symbols (numbers)
// are given in the array S, and the counts (which gives the probabilities)
// are given in C, an array of same length as S.
// x must be a number in S.
// example with symbols 0 and 1, where probabilites are P{0}=0.8, P{1}=0.2
// PutS(x,[0,1],[4,1]) and x=GetS([0,1],[4,1])
// An idea: perhaps the array S should be only in the calling function, and
// that it do not need to be passed as an argument at all.
// Hint: it may be best to to the most likely symbols (highest counts) in
// the beginning of the tables S and C.
function PutS(x,S,C)
global y Byte BitPos
global high low range ub hc lc sc K code
N=length(S);     // also length(C)
m=find(S==x);    // m is a single value, index in S, 1 <= m <= N
sc=sum(C);
lc=sc-sum(C(1:m));
hc=lc+C(m);
// disp(['PutS: lc=',int2str(lc),' hc=',int2str(hc),' sc=',int2str(sc),' m=',int2str(m)]);
EncodeSymbol;  // code the bit
return
endfunction

function x=GetS(S,C)
global y Byte BitPos
global high low range ub hc lc sc K code
range=high-low+1;
sc=sum(C);
counts=floor(( (code-low+1)*sc-1 )/range);
m=1;
lc=sc-C(1);
while (lc>counts); m=m+1; lc=lc-C(m); end;
hc=lc+C(m);
x=S(m);
// disp(['GetS: lc=',int2str(lc),' hc=',int2str(hc),' sc=',int2str(sc),' m=',int2str(m)]);
RemoveSymbol;
return
endfunction

// Aritmetic coding of a number x, 0<=x<=N, P{0}=P{1}=...=P{N}=1/(N+1)
function PutN(x,N)     // 0<=x<=N
global y Byte BitPos
global high low range ub hc lc sc K code
sc=N+1;
lc=x;
hc=x+1;
EncodeSymbol;  // code the bit
return
endfunction

function x=GetN(N)
global y Byte BitPos
global high low range ub hc lc sc K code
range=high-low+1;
sc=N+1;
x=floor(( (code-low+1)*sc-1 )/range);
hc=x+1;lc=x;
RemoveSymbol;
return
endfunction

// Aritmetic coding of a bit, probability is 0.5 for both 1 and 0
function PutABit(Bit)
global y Byte BitPos
global high low range ub hc lc sc K code
sc=2;
if Bit
    hc=1;lc=0;
else
    hc=2;lc=1;
end
EncodeSymbol;  // code the bit
return
endfunction

function Bit=GetABit
global y Byte BitPos
global high low range ub hc lc sc K code
range=high-low+1;
sc=2;
counts=floor(( (code-low+1)*sc-1 )/range);
if (1>counts)
    Bit=1;hc=1;lc=0;
else
    Bit=0;hc=2;lc=1;
end
RemoveSymbol;
return;
endfunction

// The EncodeSymbol function encode a symbol, (correspond to encode_symbol page 149)
function EncodeSymbol
global y Byte BitPos
global high low range ub hc lc sc K code
range=high-low+1;
high=low+floor(((range*hc)/sc)-1);
low=low+floor((range*lc)/sc);
while 1          // for loop on page 149
    if bitget(high,K)==bitget(low,K)
        PutBit(bitget(high,K));
        while ub > 0
            PutBit(~bitget(high,K));
            ub=ub-1;
        end
    elseif (bitget(low,K-1) & (~bitget(high,K-1)))
        ub=ub+1;
        low=bitset(low,K-1,0);
        high=bitset(high,K-1,1);
    else
        break
    end
    low=bitset(low*2,K+1,0);
    high=bitset(high*2+1,K+1,0);
end
return
endfunction

// The RemoveSymbol function removes (and fill in new) bits from
// file, y, to code
function RemoveSymbol
global y Byte BitPos
global high low range ub hc lc sc K code
range=high-low+1;
high=low+floor(((range*hc)/sc)-1);
low=low+floor((range*lc)/sc);
while 1
    if bitget(high,K)==bitget(low,K)
        // do nothing (shift bits out)
    elseif (bitget(low,K-1) & (~bitget(high,K-1)))
        code=bitset(code,K-1,~bitget(code,K-1));     // switch bit K-1
        low=bitset(low,K-1,0);
        high=bitset(high,K-1,1);
    else
        break
    end
    low=bitset(low*2,K+1,0);
    high=bitset(high*2+1,K+1,0);
    code=bitset(code*2+GetBit,K+1,0);
end
if (low > high); error('low > high'); end;
return
endfunction

// Functions to write and read a Bit
function PutBit(Bit)
global y Byte BitPos
BitPos=BitPos-1;
if (~BitPos); Byte=Byte+1; BitPos=8; end;
y(Byte) = bitset(y(Byte),BitPos,Bit);
return
endfunction

function Bit=GetBit
global y Byte BitPos
BitPos=BitPos-1;
if (~BitPos); Byte=Byte+1; BitPos=8; end;
Bit=bitget(y(Byte),BitPos);
return
endfunction

function b=BitEst(N,N1);
if (N1>(N/2)); N1=N-N1; end;
N0=N-N1;
if (N>1000)
    b=(N+3/2)*log2(N)-(N0+1/2)*log2(N0)-(N1+1/2)*log2(N1)-1.3256;
elseif (N1>20)
    b=(N+3/2)*log2(N)-(N0+1/2)*log2(N0)-(N1+1/2)*log2(N1)-0.020984*log2(log2(N))-1.25708;
else
    b=log2(N+1)+sum(log2(N-(0:(N1-1))))-sum(log2(N1-(0:(N1-1))));
end
return
endfunction