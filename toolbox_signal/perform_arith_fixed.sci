function y = perform_arith_fixed(x, h, n)

// perform_arith_fixed - arithmetic coding
// 
// Coding:
//   y = perform_arith_fixed(x, h, n);
// Decoding:
//   x = perform_arith_fixed(y, h);
//
//   Copyright (c) 2008 Gabriel Peyre

if argn(2)==2
    direction=1;
elseif argn(2)==3
    direction=-1;
else
    error('Wrong number of arumgnets');
end

h = round(h*100000); 

if direction==1
    y = arithenco(x,h);
else
    y = arithdeco(x, h, n);
end
endfunction

function code = arithenco(seq, counts)
//ARITHENCO Encode a sequence of symbols using arithmetic coding.
//   CODE = ARITHENCO(SEQ, COUNTS) generates binary arithmetic code 
//   corresponding to the sequence of symbols specified in the vector SEQ. 
//   The vector COUNTS contains the symbol counts (the number of times each
//   symbol of the source's alphabet occurs in a test data set) and represents
//   the source's statistics.
//
//   Example: 
//     Consider a source whose alphabet is {x, y, z}. A 177-symbol test data 
//     set from the source contains 29 x's, 48 y's and 100 z's. To encode the 
//     sequence yzxzz, use these commands:
//
//       seq = [2 3 1 3 3];
//       counts = [29 48 100];
//       code = arithenco(seq, counts) 
//                   
//   See also ARITHDECO.

//   Copyright 1996-2002 The MathWorks, Inc.
//   $Revision: 1.3 $ $Date: 2002/06/17 12:22:10 $

//   References:
//         [1] Sayood, K., Introduction to Data Compression, 
//         Morgan Kaufmann, 2000, Chapter 4, Section 4.4.3.
         


// Check the incoming orientation and adjust if necessary
[row_s, col_s] = size(seq);
if (row_s > 1),
    seq = seq';
end

[row_c, col_c] = size(counts);
if (row_c > 1),
    counts = counts';
end

// Compute the cumulative counts vector from the counts 
cum_counts = [0, cumsum(counts)];

// Compute the Word Length required.
total_count = cum_counts($);
N = ceil(log2(total_count)) + 2;

// Initialize the lower and upper bounds.
dec_low = 0;
dec_up = 2^N-1;
E3_count = 0;

// Obtain an over estimate for the length of CODE and initialize CODE
code_len = length(seq) * ( ceil(log2(length(counts))) + 2 ) + N;
code = zeros(1, code_len);
code_index = 1;

// Loop for each symbol in SEQ
for k = 1:length(seq)

    symbol = seq(k);
    // Compute the new lower bound
    dec_low_new = dec_low + floor( (dec_up-dec_low+1)*cum_counts(symbol+1-1)/total_count );

    // Compute the new upper bound
    dec_up = dec_low + floor( (dec_up-dec_low+1)*cum_counts(symbol+1)/total_count )-1;

    // Update the lower bound
    dec_low = dec_low_new;
    
    // Check for E1, E2 or E3 conditions and keep looping as long as they occur.
    while( isequal(bitget(dec_low, N), bitget(dec_up, N)) | ..
        (isequal(bitget(dec_low, N-1), 1) & isequal(bitget(dec_up, N-1), 0) ) ),
        
        // If it is an E1 or E2 condition,
        if isequal(bitget(dec_low, N), bitget(dec_up, N)),

            // Get the MSB
            b = bitget(dec_low, N);
            code(code_index) = b;
            code_index = code_index + 1;
        
            // Left shifts
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up = bitshift(dec_up, 1) + 1;
            
            // Check if E3_count is non-zero and transmit appropriate bits
            if (E3_count > 0),
                // Have to transmit complement of b, E3_count times.
                code(code_index:code_index+E3_count-1) = bitcmp(b, 1).*ones(1, E3_count);
                code_index = code_index + E3_count;
                E3_count = 0;
            end

            // Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up, N+1, 0);
            
        // Else if it is an E3 condition    
        elseif ( (isequal(bitget(dec_low, N-1), 1) & ...
            isequal(bitget(dec_up, N-1), 0) ) ),
            
            // Left shifts
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up  = bitshift(dec_up, 1) + 1;

            // Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up, N+1, 0);
            
            // Complement the new MSB of dec_low and dec_up
            dec_low = bitxor(dec_low, 2^(N-1) );
            dec_up  = bitxor(dec_up, 2^(N-1) );
            
            // Increment E3_count to keep track of number of times E3 condition is hit.
            E3_count = E3_count+1;
        end
    end
end
 
// Terminate encoding
bin_low = de2bi(dec_low, N, 'left-msb');
if E3_count==0,
    // Just transmit the final value of the lower bound bin_low       
    code(code_index:code_index + N - 1) = bin_low;
    code_index = code_index + N;
else
   // Transmit the MSB of bin_low. 
   b = bin_low(1);
   code(code_index) = b;
   code_index = code_index + 1;
   
   // Then transmit complement of b (MSB of bin_low), E3_count times. 
   code(code_index:code_index+E3_count-1) = bitcmp(b, 1).*ones(1, E3_count);
   code_index = code_index + E3_count;

   // Then transmit the remaining bits of bin_low
   code(code_index:code_index+N-2) = bin_low(2:N);
   code_index = code_index + N - 1;
end          

// Output only the filled values
code = code(1:code_index-1);

// Set the same output orientation as seq
if (row_s > 1)
    code = code.';
end



endfunction

function dseq = arithdeco(code, counts, len)
//ARITHDECO Decode binary code using arithmetic decoding.
//   DSEQ = ARITHDECO(CODE, COUNTS, LEN) decodes the binary arithmetic code
//   in the vector CODE (generated using ARITHENCO) to the corresponding
//   sequence of symbols. The vector COUNTS contains the symbol counts (the
//   number of times each symbol of the source's alphabet occurs in a test
//   data set) and represents the source's statistics. LEN is the number of
//   symbols to be decoded. 
//   
//   Example: 
//     Consider a source whose alphabet is {x, y, z}. A 177-symbol test data 
//     set from the source contains 29 x's, 48 y's and 100 z's. To encode the 
//     sequence yzxzz, use these commands:
//
//       seq = [2 3 1 3 3];
//       counts = [29 48 100];
//       code = arithenco(seq, counts)   
//            
//     To decode this code (and recover the sequence of  
//     symbols it represents) use this command:
//            
//       dseq = arithdeco(code, counts, 5)
//            
//   See also ARITHENCO.

//   Copyright 1996-2002 The MathWorks, Inc.
//   $Revision: 1.2 $ $Date: 2002/04/14 20:12:32 $

//   References:
//         [1] Sayood, K., Introduction to Data Compression, 
//         Morgan Kaufmann, 2000, Chapter 4, Section 4.4.3.

// Check the incoming orientation and adjust if necessary
[row_cd, col_cd] = size(code);
if (row_cd > 1),
    code = code';
end

[row_c, col_c] = size(counts);
if (row_c > 1),
    counts = counts';
end

// Compute the cumulative counts vector from the counts vector
cum_counts = [0, cumsum(counts)];

// Compute the Word Length (N) required.
total_count = cum_counts($);
N = ceil(log2(total_count)) + 2;

// Initialize the lower and upper bounds.
dec_low = 0;
dec_up = 2^N-1;

// Read the first N number of bits into a temporary tag bin_tag
bin_tag = code(1:N);
dec_tag = bi2de(bin_tag, 'left-msb');

// Initialize DSEQ
dseq = zeros(1,len);
dseq_index = 1;

k=N;
ptr = 0;

// This loop runs untill all the symbols are decoded into DSEQ
while (dseq_index <= len)
    
    // Compute dec_tag_new
    dec_tag_new =floor( ((dec_tag-dec_low+1)*total_count-1)/(dec_up-dec_low+1) );
    
    // Decode a symbol based on dec_tag_new
    ptr = pick(cum_counts, dec_tag_new);
    
    // Update DSEQ by adding the decoded symbol
    dseq(dseq_index) = ptr;
    dseq_index = dseq_index + 1;
    
    // Compute the new lower bound
    dec_low_new = dec_low + floor( (dec_up-dec_low+1)*cum_counts(ptr-1+1)/total_count );
    
    // Compute the new upper bound
    dec_up = dec_low + floor( (dec_up-dec_low+1)*cum_counts(ptr+1)/total_count )-1;
    
    // Update the lower bound
    dec_low = dec_low_new;
    
    // Check for E1, E2 or E3 conditions and keep looping as long as they occur.
     while ( isequal(bitget(dec_low, N), bitget(dec_up, N)) | ...
        ( isequal(bitget(dec_low, N-1), 1) & isequal(bitget(dec_up, N-1), 0) ) ),
        
        // Break out if we have finished working with all the bits in CODE
        if ( k==length(code) ), break, end;
        k = k + 1;

        // If it is an E1 or E2 condition, do
        if isequal(bitget(dec_low, N), bitget(dec_up, N)),

            // Left shifts and update
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up  = bitshift(dec_up,  1) + 1;

            // Left shift and read in code
            dec_tag = bitshift(dec_tag, 1) + code(k);

            // Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up,  N+1, 0);
            dec_tag = bitset(dec_tag, N+1, 0);
        
        // Else if it is an E3 condition        
        elseif ( isequal(bitget(dec_low, N-1), 1) & ...
            isequal(bitget(dec_up, N-1), 0) ),

            // Left shifts and update
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up  = bitshift(dec_up,  1) + 1;

            // Left shift and read in code
            dec_tag = bitshift(dec_tag, 1) + code(k);
            
            // Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up,  N+1, 0);
            dec_tag = bitset(dec_tag, N+1, 0);
            
            // Complement the new MSB of dec_low, dec_up and dec_tag
            dec_low = bitxor(dec_low, 2^(N-1) );
            dec_up  = bitxor(dec_up,  2^(N-1) );
            dec_tag = bitxor(dec_tag, 2^(N-1) );

        end
    end // end while
end // end while length(dseq)

// Set the same output orientation as code
if (row_cd > 1)
    dseq = dseq.';
end


endfunction

//-------------------------------------------------------------
function [ptr] = pick(cum_counts, value);
// This internal function is used to find where value is positioned

// Check for this case and quickly exit
if value == cum_counts($)
    ptr = length(cum_counts)-1;
    return
end

c = find(cum_counts <= value);
ptr = c($);


endfunction


function y = bi2de(x, s)

x = x(:)';
n = length(x);
s = (0:n-1);
s = s($:-1:1);
y = sum( x.* 2.^s );

endfunction


// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) ???? - INRIA - Farid BELAHCENE
// Copyright (C) 2006 - INRIA - Pierre MARECHAL
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// =============================================================================
//
// BIN2dec function
// Given str a binary string, this function returns the decimal number whose the
// binary representation is given by str
//
// -Input :
//    str : a string (or a vector/matrix of strings), containing only characters
//         '1' and '0'
// -Output :
//    y : a scalar/vector/matrix
//
// F.Belahcene

// check the type of input argument

// 2006-06-26 : Modified by Pierre MARECHAL
// Check length of given string ( must be 47 bits or less )
// =============================================================================

function y=bin2dec(str)

	if type(str)<>10
		error(msprintf(gettext("%s: Wrong type for input argument #%d: Matrix of strings expected.\n"),"bin2dec",1));
	end
	
	// delete all spaces included in the str
	str = strsubst(str,' ','');
	
	// check the str characters are only '0' or '1', and convert the binary str to corresponing decimal number
	for i=1:prod(size(str))
		
		ind1=strindex(str(i),'1')
		ind0=strindex(str(i),'0')
		
		if length(str(i)) <> sum([prod(size(ind0)) prod(size(ind1))]) then
			error(msprintf(gettext("%s: Wrong value for input argument #%d: Matrix of strings made of zeros and ones expected.\n"),"bin2dec",1));
		end
		
		if length(str(i)) > 54 then
			error(msprintf(gettext("%s: Wrong size for input argument #%d: Must be less than %d characters.\n"),"bin2dec",1,54));
		end
		
		if ~isempty(ind1)
			ind1   = length(str(i))-ind1($:-1:1);
			y($+1) = sum(2^ind1);
		elseif ~isempty(ind0)
			y($+1) = 0;
		else
			y($+1) = [];
		end
	end
	
	y=matrix(y,size(str));
	
endfunction


function gettext(x)
endfunction


function X = bitshift (A,k,n)
X = A*2^k;
endfunction


function z=bitxor(x,y)
// Copyright INRIA
// BITXOR function
// Given x,y two positives integers this function returns the decimal number whose the binary form is the XOR of the binary representations of x and y

// -Inputs : 
//  x, y :  scalars/vectors/matices/hypermatices of positives integers, x and y must have the same size
// -Output :
//  z : scalar/vector/matrix/hypermatice
//
// F.Belahcene

if typeof(x)<>typeof(y)
	error('Error');
end
if typeof(x)<>"uint8" & typeof(x)<>"uint16" & typeof(x)<>"uint32" & typeof(x)<>"constant"
	error('Error');
end

if size(x)<>size(y)
	error('Error');
elseif isempty(x) & isempty(x)
	z=[]
	return
end

if (type(x)==1 & (x-floor(x)<>0 | x<0)) | (type(x)==8 & x<0) | (type(x)==17 & (type(x.entries<>1) | type(x.entries<>8)) & find(x.entries>0)<>[])  | (type(x)<>1 & type(x)<>8 & type(x)<>17)
	error('Error');
end

if (type(y)==1 & (y-floor(y)<>0 | y<0)) | (type(y)==8 & y<0) | (type(y)==17 & (type(y.entries<>1) | type(y.entries<>8)) & find(y.entries>0)<>[]) | (type(y)<>1 & type(y)<>8 & type(y)<>17)
	error('Error');
end

for i=1:prod(size(x))
	zbin=dec2bin([x(i);y(i)])
	zand=strcat((string(sum(asciimat(zbin)-48,1))))
	zand=strsubst(zand,'2','0')
	z(i)=bin2dec(zand)
end
z=matrix(z,size(x))

if typeof(x)=="uint8"
	z=uint8(z)
elseif typeof(x)=="uint16"
	z=uint16(z)
elseif typeof(x)=="uint32"
	z=uint32(z)
end

endfunction


function y = de2bi(x,N,str)

y = str2code(dec2bin(x,N))';

endfunction


// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) INRIA - Farid BELAHCENE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// =============================================================================
//
// Author : F.Belahcene
// DEC2BIN function
//
// Given x, a positive scalar/vector/matix of reals, this function returns
// a column vector of strings. Each string is the binary representation
// of the input argument components (i.e y(i) is the binary representation
// of x(i))
//
// -Inputs :
//    x : a  scalar/vector/matix of positives reals
//    n : an integer
// -Output :
//    y : a vector of strings (positives)
//
// =============================================================================

function y=dec2bin(x,n)
	
	rhs = argn(2);
	
	// check the number of input arguments
	if rhs<1 or rhs>2 then
		error();
	end
	
	// check type and size of the input arguments
	if or(type(x)<>8) & (or(type(x)<>1) | or(x<0)) then
		error();
	end
	
	if rhs==2 & ((type(n)<>8 & (type(n)<>1 | n<0)) | prod(size(n))<>1) then
		error();
	end
	
	// empty matrix
	if x==[]
		y=string([]);
		return;
	end
	
	[nr,nc] = size(x);
	x=x(:);
	
	// input argument is a scalar/vector/matrix of zeros
	
	if and(x==0)
		if rhs==2
			y = strcat(string(zeros(1:n))) + emptystr(nr,nc);
		else
			y = "0" + emptystr(nr,nc);
		end
		return
	end
	
	// for x=25, pow=[4 3 0], because x=2^4+2^3+2^0
	// for x=23, pow=[4 2 1 0] because x=2^4+2^2+2^1+2^0
	// for x=[25 23]
	// pow=[4 3 0 -1
	//      4 2 1 0];
	
	while find(x>0)<>[]
		pow(x>0,$+1) = floor(log2(double(x(x>0))));
		pow(x<=0,$)  = -1;
		x(x>0)       = floor(x(x>0)-2^pow(x>0,$));
	end
	
	pow   = pow+1;
	ytemp = zeros(size(pow,1),size(pow,2));
	
	for i=1:size(ytemp,1)
		ind          = pow(i,pow(i,:)>=1);
		ytemp(i,ind) = 1;
	end
	
	if rhs==2
		for i=1:size(ytemp,1)
			y(i)=strcat(string([zeros(1,round(n-size(ytemp,2))) ytemp(i,size(ytemp,2):-1:1)]));
		end
	else
		for i=1:size(ytemp,1)
			y(i)=strcat(string(ytemp(i,size(ytemp,2):-1:1)));
		end
	end
	
	y = matrix(y,nr,nc);
	
endfunction



// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) ???? - INRIA - Farid BELAHCENE
// Copyright (C) 2008 - INRIA - Pierre MARECHAL
// 
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function y = bitcmp(x,n)
	
	// BITCMP function
	//
	// Given an unsigned integer x, this function returns the unsigned integer
	// which is the integer corresponding to the complementary of the binary
	// form of x
	//
	// If the bits number of the x binary representation is less than the
	// bitmax number (8,16 or 32) then the bits '1' are added to the
	// complementary in order to have bitmax number (8, 16 or 32) for the
	// complementary
	//
	// for example for the type uint8 (bitmax=8), the complementary of '1101' is not '0010' but '11110010'
	// The integer n sets the bits max number
	// -Inputs :
	//  x : an unsigned integer
	//  n : a positive integer between 1 and the bitmax of the x type
	//
	// -Output :
	//  y : an unsigned integer
	//
	// P. Marechal, 5 Feb 2008
	//   - Add argument check
	
	// Check input arguments
	// =========================================================================
	
	rhs = argn(2);
	
	// check number input argument
	
	rhs = argn(2);
	if (type(x) == 1) & (rhs == 1) then
        error('Problem');
	elseif rhs == 0 then
        error('Problem');
	end
	
	// check type
	
	if    (type(x)==1  & (x-floor(x)<>0 | x<0)) ..
		| (type(x)==8  & (inttype(x)<10)) ..
		| (type(x)<>1  & type(x)<>8) then
		
        error('Problem');
	end
	
	if  (rhs == 2) & ( ..
			(type(n)==1  & (n-floor(n)<>0 | x<0)) ..
			| (type(n)==8  & (inttype(n)<10)) ..
			| (type(n)<>1  & type(n)<>8) ..
			| (size(n,"*")<>1) ) then
        error('Problem');
	end
	
	// check n value
	
	if rhs>1 then
		
		select inttype(x)
			case 0  then nmax = 52;
			case 11 then nmax = 8;
			case 12 then nmax = 16;
			case 14 then nmax = 32;
		end
		
		if (n>nmax) | (n<1) then
        error('Problem');
		end
		
	else
		n = nmax;
	end
	
	// Algorithm
	// =========================================================================
	
	// empty matrix shortcut
	
	if isempty(x) then
		y = [];
		return;
	end
	
	// unit8, uint16 and uint32 shortcut
	
	if type(x)==8 then
		y = ~x;
		if rhs > 1 then
			select inttype(x)
				case 11 then y = y & uint8(  2^n - 1);
				case 12 then y = y & uint16( 2^n - 1);
				case 14 then y = y & uint32( 2^n - 1);
			end
		end
		return;
	end
	
	n = ones(x)*n;
	
	if type(x) == 1 then
		
		a     = 2^32;
		
		y_LSB = uint32( x - double(uint32(x/a)) * a ); // LSB Less Significant Bits
		y_MSB = uint32( x/a );                         // MSB Most Significant Bits
		
		y_LSB = ~y_LSB;
		y_MSB = ~y_MSB;
		
		if n <= 32 then
			y_LSB = y_LSB & uint32( 2^n - 1);
			y_MSB = uint32(0);
		else
			y_MSB = y_MSB & uint32( 2^(n-32) - 1);
		end
		
		y = double( a * y_MSB + y_LSB );
	end
	
endfunction

// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) ???? - INRIA - Farid BELAHCENE
// Copyright (C) 2008 - INRIA - Pierre MARECHAL
// 
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function y = bitget(x,pos)
	
	// BITGET function
	// Given an unsigned integer x, this function returns an unsigned integer
	// (0 or 1) which is the bit number 'pos' from the representation binary of x
	// if x=uint8(19) and pos=2 then bitget returns the 2th bit of the binary form of 19 ('10011') which is 1
	// -Inputs :
	//  x : an unsigned integer
	// pos : a positive integer between 1 and the bitmax of the x type
	//
	// -Output :
	//  y : an unsigned integer
	//
	// F.Belahcene
	
	// P. Marechal, 5 Feb 2008
	//   - Add argument check
	
	// Check input arguments
	// =========================================================================
	
	// check number input argument
	
	rhs = argn(2);
	
	if rhs <> 2 then
        error('Problem');
	end
	
	// case empty matrix
	
	if isempty(x)
		if ~isempty(pos) & prod(size(pos))<>1
        error('Problem');
		else
			y=[]
			return
		end
	end
	
	// check size
	
	if (size(x,"*")>1) & (size(pos,"*")>1) & (or(size(x)<>size(pos))) then
        error('Problem');
	end
	
	// check type
	
	if    (type(x)==1  & (x-floor(x)<>0 | x<0)) ..
		| (type(x)==8  & (inttype(x)<10)) ..
		| (type(x)<>1  & type(x)<>8) then
		
        error('Problem');
	end
	
	if    (type(pos)==1  & (pos-floor(pos)<>0 | pos<0)) ..
		| (type(pos)==8  & (inttype(pos)<10)) ..
		| (type(pos)<>1  & type(pos)<>8) then
		
        error('Problem');
	end
	
	// check pos value
	
	select inttype(x)
		case 0  then posmax = 52;
		case 11 then posmax = 8;
		case 12 then posmax = 16;
		case 14 then posmax = 32;
	end
	
	if (pos>posmax) | (pos<1) then
        error('Problem');
	end
	
	// Algorithm
	// =========================================================================
	
	if size(pos,"*") == 1;
		pos  = ones(x)  * pos;
	end
	
	if size(x,"*") == 1;
		x    = ones(pos) * x;
	end
	
	if type(x)==8 then
		
		select inttype(x)
			
			case 11 then
				mask = uint8(2^(pos-1));
				y    = uint8(1 * ((x & mask) > 0));
				return;
			
			case 12 then
				mask = uint16(2^(pos-1));
				y    = uint16(1 * ((x & mask) > 0));
				return;
			
			case 14 then
				mask = uint32(2^(pos-1));
				y    = uint32(1 * ((x & mask) > 0));
				return;
		end
		
	else
		
		// type == 1
		
		a     = 2^32;
		mask  = uint32(zeros(pos));
		ytemp = uint32(zeros(pos));
		
		if or(pos<=32) then
			mask( pos<=32 )  = uint32( 2^(pos(pos<=32) -1 ));
			ytemp( pos<=32 ) = uint32( x(pos<=32) - double(uint32(x(pos<=32)/a)) * a ); // permet de r?cup?rer les 32 bits de poid faible
		end
		
		if or(pos>32) then
			mask( pos>32  )     = uint32( 2^(pos(pos>32) -32 -1 ));
			ytemp( pos> 32 )    = uint32( x(pos> 32)/a); // permet de r?cup?rer les 32 bits de poid fort
		end
		
		y = 1 * ((ytemp & mask) > 0);
		
	end
	
endfunction


// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2008 - INRIA - Pierre MARECHAL
// 
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

function y = bitset(x,pos,v)
	
	// INRIA 2008 - Pierre MARECHAL
	//
	// BITSET function
	// Set bit at specified position
	
	// Check input arguments
	// =========================================================================
	
	// check number input argument
	
	rhs = argn(2);
	
	if rhs < 2 then
		error('Error');
	end
	
	// case empty matrix
	
	if isempty(x)
		if ~isempty(pos) & prod(size(pos))<>1
			error('Error');
		else
			y=[]
			return
		end
	end
	
	// check size
	
	if (size(x,"*")>1) & (size(pos,"*")>1) & (or(size(x)<>size(pos))) then
		error('Error');
	end
	
	// check type
	
	if    (type(x)==1  & (x-floor(x)<>0 | x<0)) ..
		| (type(x)==8  & (inttype(x)<10)) ..
		| (type(x)<>1  & type(x)<>8) then
		error('Error');
	end
	
	if    (type(pos)==1  & (pos-floor(pos)<>0 | pos<0)) ..
		| (type(pos)==8  & (inttype(pos)<10)) ..
		| (type(pos)<>1  & type(pos)<>8) then
		error('Error');
	end
	
	// check pos value
	
	select inttype(x)
		case 0  then posmax = 52;
		case 11 then posmax = 8;
		case 12 then posmax = 16;
		case 14 then posmax = 32;
	end
	
	if (pos>posmax) | (pos<1) then
        error('Error');
	end
	
	// check v value
	
	if rhs == 3 & ..
		( ((type(v)<>1) & (type(v)<>8)) ..
		| ((type(x)==8) & (inttype(x)<10)) ..
		| ((v<>0) & (v<>1)) ) then
		error('Error');
	end
	
	// Algorithm
	// =========================================================================
	
	if size(pos,"*") == 1;
		pos  = ones(x)*pos;
	end
	
	if size(x,"*") == 1;
		x    = ones(pos)*x;
	end
	
	if rhs<3 then
		v = 1;
	end
	
	if type(x)==8 then
		
		select inttype(x)
			case 11 then mask = uint8(2^(pos-1));
			case 12 then mask = uint16(2^(pos-1));
			case 14 then mask = uint32(2^(pos-1));
		end
		
		if v==0 then
			y = x & (~mask);
		else
			y = x | mask;
		end
		
		return;
		
	else
		
		// type == 1
		
		a     = 2^32;
		mask  = uint32(zeros(pos));
		
		y_MSB  = uint32(zeros(pos));
		y_LSB  = uint32(zeros(pos));
		
		y_LSB = uint32( x - double(uint32(x/a)) * a ); // LSB Less Significant Bits
		y_MSB = uint32( x/a );                         // MSB Most Significant Bits
		
		if or(pos<=32) then
			mask(  pos<=32 ) = uint32( 2^(pos(pos<=32) -1 ));
			if v==0 then
				y_LSB( pos<=32 ) = y_LSB(pos<=32) & (~mask(pos<=32));
			else
				y_LSB( pos<=32 ) = y_LSB(pos<=32) | mask(pos<=32);
			end
		end
		
		if or(pos>32) then
			mask(  pos>32 ) = uint32( 2^(pos(pos>32) -32 -1 ));
			if v==0 then
				y_MSB( pos>32 ) = y_MSB(pos>32) & (~ mask(pos>32));
			else
				y_MSB( pos>32 ) = y_MSB(pos>32) | mask(pos>32);
			end
		end
		
		y = double( a * y_MSB + y_LSB );
		
	end
	
endfunction

