function y = perform_arith_fixed(x, h, n)

% perform_arith_fixed - arithmetic coding
% 
% Coding:
%   y = perform_arith_fixed(x, h, n);
% Decoding:
%   x = perform_arith_fixed(y, h);
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin==2
    dir=1;
elseif nargin==3
    dir=-1;
else
    error('Wrong number of arumgnets');
end

h = round(h*100000); 

if dir==1
    y = arithenco(x,h);
else
    y = arithdeco(x, h, n);
end



function code = arithenco(seq, counts)
%ARITHENCO Encode a sequence of symbols using arithmetic coding.
%   CODE = ARITHENCO(SEQ, COUNTS) generates binary arithmetic code 
%   corresponding to the sequence of symbols specified in the vector SEQ. 
%   The vector COUNTS contains the symbol counts (the number of times each
%   symbol of the source's alphabet occurs in a test data set) and represents
%   the source's statistics.
%
%   Example: 
%     Consider a source whose alphabet is {x, y, z}. A 177-symbol test data 
%     set from the source contains 29 x's, 48 y's and 100 z's. To encode the 
%     sequence yzxzz, use these commands:
%
%       seq = [2 3 1 3 3];
%       counts = [29 48 100];
%       code = arithenco(seq, counts) 
%                   
%   See also ARITHDECO.

%   Copyright 1996-2002 The MathWorks, Inc.
%   $Revision: 1.3 $ $Date: 2002/06/17 12:22:10 $

%   References:
%         [1] Sayood, K., Introduction to Data Compression, 
%         Morgan Kaufmann, 2000, Chapter 4, Section 4.4.3.
         

% Check the incoming orientation and adjust if necessary
[row_s, col_s] = size(seq);
if (row_s > 1),
    seq = seq.';
end

[row_c, col_c] = size(counts);
if (row_c > 1),
    counts = counts.';
end

% Compute the cumulative counts vector from the counts 
cum_counts = [0, cumsum(counts)];

% Compute the Word Length required.
total_count = cum_counts(end);
N = ceil(log2(total_count)) + 2;

% Initialize the lower and upper bounds.
dec_low = 0;
dec_up = 2^N-1;
E3_count = 0;

% Obtain an over estimate for the length of CODE and initialize CODE
code_len = length(seq) * ( ceil(log2(length(counts))) + 2 ) + N;
code = zeros(1, code_len);
code_index = 1;

% Loop for each symbol in SEQ
for k = 1:length(seq)

    symbol = seq(k);
    % Compute the new lower bound
    dec_low_new = dec_low + floor( (dec_up-dec_low+1)*cum_counts(symbol+1-1)/total_count );

    % Compute the new upper bound
    dec_up = dec_low + floor( (dec_up-dec_low+1)*cum_counts(symbol+1)/total_count )-1;

    % Update the lower bound
    dec_low = dec_low_new;
    
    % Check for E1, E2 or E3 conditions and keep looping as long as they occur.
    while( isequal(bitget(dec_low, N), bitget(dec_up, N)) || ...
        (isequal(bitget(dec_low, N-1), 1) && isequal(bitget(dec_up, N-1), 0) ) ),
        
        % If it is an E1 or E2 condition,
        if isequal(bitget(dec_low, N), bitget(dec_up, N)),

            % Get the MSB
            b = bitget(dec_low, N);
            code(code_index) = b;
            code_index = code_index + 1;
        
            % Left shifts
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up = bitshift(dec_up, 1) + 1;
            
            % Check if E3_count is non-zero and transmit appropriate bits
            if (E3_count > 0),
                % Have to transmit complement of b, E3_count times.
                code(code_index:code_index+E3_count-1) = bitcmp(b, 1).*ones(1, E3_count);
                code_index = code_index + E3_count;
                E3_count = 0;
            end

            % Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up, N+1, 0);
            
        % Else if it is an E3 condition    
        elseif ( (isequal(bitget(dec_low, N-1), 1) && ...
            isequal(bitget(dec_up, N-1), 0) ) ),
            
            % Left shifts
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up  = bitshift(dec_up, 1) + 1;

            % Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up, N+1, 0);
            
            % Complement the new MSB of dec_low and dec_up
            dec_low = bitxor(dec_low, 2^(N-1) );
            dec_up  = bitxor(dec_up, 2^(N-1) );
            
            % Increment E3_count to keep track of number of times E3 condition is hit.
            E3_count = E3_count+1;
        end
    end
end
 
% Terminate encoding
bin_low = de2bi(dec_low, N, 'left-msb');
if E3_count==0,
    % Just transmit the final value of the lower bound bin_low       
    code(code_index:code_index + N - 1) = bin_low;
    code_index = code_index + N;
else
   % Transmit the MSB of bin_low. 
   b = bin_low(1);
   code(code_index) = b;
   code_index = code_index + 1;
   
   % Then transmit complement of b (MSB of bin_low), E3_count times. 
   code(code_index:code_index+E3_count-1) = bitcmp(b, 1).*ones(1, E3_count);
   code_index = code_index + E3_count;

   % Then transmit the remaining bits of bin_low
   code(code_index:code_index+N-2) = bin_low(2:N);
   code_index = code_index + N - 1;
end          

% Output only the filled values
code = code(1:code_index-1);

% Set the same output orientation as seq
if (row_s > 1)
    code = code.';
end

%----------------------------------------------------------
function eStr = errorchk(seq, counts)
% Function for validating the input parameters. 

eStr.ecode = 0;
eStr.emsg = '';

% Check to make sure a vector has been entered as input and not a matrix
if (length(find(size(seq)==1)) ~= 1)
    eStr.emsg = ['The symbol sequence parameter must be a vector of positive ',...
                 'finite integers.'];
    eStr.ecode = 1; return;    
end

% Check to make sure a character array is not specified for SEQ
if ischar(seq)
    eStr.emsg = ['The symbol sequence parameter must be a vector of positive ',...
                 'finite integers.'];
    eStr.ecode = 1; return;    
end

% Check to make sure that finite positive integer values (non-complex) are 
% entered for SEQ
if ~all(seq > 0) || ~all(isfinite(seq)) || ~isequal(seq, round(seq)) || ...
    ~isreal(seq)
    eStr.emsg = ['The symbol sequence parameter must be a vector of positive ',...
                 'finite integers.'];
    eStr.ecode = 1; return;    
end

if length(find(size(counts)==1)) ~= 1
    eStr.emsg = ['The symbol counts parameter must be a vector of positive ',...
                 'finite integers.'];
    eStr.ecode = 1; return;    
end

% Check to make sure that finite positive integer values (non-complex) are
% entered for COUNTS
if ~all(counts > 0) || ~all(isfinite(counts)) || ~isequal(counts, round(counts)) || ...
     ~isreal(counts)
    eStr.emsg = ['The symbol counts parameter must be a vector of positive ',...
                 'finite integers.'];
    eStr.ecode = 1; return;    
end

% Check to ensure that the maximum value in the SEQ vector is less than or equal
% to the length of the counts vector COUNTS.
if max(seq) > length(counts)
    eStr.emsg = ['The symbol sequence parameter can take values only between',...
                 ' 1 and the length of the symbol counts parameter.'];
    eStr.ecode = 1; return;    
end

% [EOF]



function dseq = arithdeco(code, counts, len)
%ARITHDECO Decode binary code using arithmetic decoding.
%   DSEQ = ARITHDECO(CODE, COUNTS, LEN) decodes the binary arithmetic code
%   in the vector CODE (generated using ARITHENCO) to the corresponding
%   sequence of symbols. The vector COUNTS contains the symbol counts (the
%   number of times each symbol of the source's alphabet occurs in a test
%   data set) and represents the source's statistics. LEN is the number of
%   symbols to be decoded. 
%   
%   Example: 
%     Consider a source whose alphabet is {x, y, z}. A 177-symbol test data 
%     set from the source contains 29 x's, 48 y's and 100 z's. To encode the 
%     sequence yzxzz, use these commands:
%
%       seq = [2 3 1 3 3];
%       counts = [29 48 100];
%       code = arithenco(seq, counts)   
%            
%     To decode this code (and recover the sequence of  
%     symbols it represents) use this command:
%            
%       dseq = arithdeco(code, counts, 5)
%            
%   See also ARITHENCO.

%   Copyright 1996-2002 The MathWorks, Inc.
%   $Revision: 1.2 $ $Date: 2002/04/14 20:12:32 $

%   References:
%         [1] Sayood, K., Introduction to Data Compression, 
%         Morgan Kaufmann, 2000, Chapter 4, Section 4.4.3.


% Check the incoming orientation and adjust if necessary
[row_cd, col_cd] = size(code);
if (row_cd > 1),
    code = code.';
end

[row_c, col_c] = size(counts);
if (row_c > 1),
    counts = counts.';
end

% Compute the cumulative counts vector from the counts vector
cum_counts = [0, cumsum(counts)];

% Compute the Word Length (N) required.
total_count = cum_counts(end);
N = ceil(log2(total_count)) + 2;

% Initialize the lower and upper bounds.
dec_low = 0;
dec_up = 2^N-1;

% Read the first N number of bits into a temporary tag bin_tag
bin_tag = code(1:N);
dec_tag = bi2de(bin_tag, 'left-msb');

% Initialize DSEQ
dseq = zeros(1,len);
dseq_index = 1;

k=N;
ptr = 0;

% This loop runs untill all the symbols are decoded into DSEQ
while (dseq_index <= len)
    
    % Compute dec_tag_new
    dec_tag_new =floor( ((dec_tag-dec_low+1)*total_count-1)/(dec_up-dec_low+1) );
    
    % Decode a symbol based on dec_tag_new
    ptr = pick(cum_counts, dec_tag_new);
    
    % Update DSEQ by adding the decoded symbol
    dseq(dseq_index) = ptr;
    dseq_index = dseq_index + 1;
    
    % Compute the new lower bound
    dec_low_new = dec_low + floor( (dec_up-dec_low+1)*cum_counts(ptr-1+1)/total_count );
    
    % Compute the new upper bound
    dec_up = dec_low + floor( (dec_up-dec_low+1)*cum_counts(ptr+1)/total_count )-1;
    
    % Update the lower bound
    dec_low = dec_low_new;
    
    % Check for E1, E2 or E3 conditions and keep looping as long as they occur.
     while ( isequal(bitget(dec_low, N), bitget(dec_up, N)) | ...
        ( isequal(bitget(dec_low, N-1), 1) & isequal(bitget(dec_up, N-1), 0) ) ),
        
        % Break out if we have finished working with all the bits in CODE
        if ( k==length(code) ), break, end;
        k = k + 1;

        % If it is an E1 or E2 condition, do
        if isequal(bitget(dec_low, N), bitget(dec_up, N)),

            % Left shifts and update
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up  = bitshift(dec_up,  1) + 1;

            % Left shift and read in code
            dec_tag = bitshift(dec_tag, 1) + code(k);

            % Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up,  N+1, 0);
            dec_tag = bitset(dec_tag, N+1, 0);
        
        % Else if it is an E3 condition        
        elseif ( isequal(bitget(dec_low, N-1), 1) & ...
            isequal(bitget(dec_up, N-1), 0) ),

            % Left shifts and update
            dec_low = bitshift(dec_low, 1) + 0;
            dec_up  = bitshift(dec_up,  1) + 1;

            % Left shift and read in code
            dec_tag = bitshift(dec_tag, 1) + code(k);
            
            % Reduce to N for next loop
            dec_low = bitset(dec_low, N+1, 0);
            dec_up  = bitset(dec_up,  N+1, 0);
            dec_tag = bitset(dec_tag, N+1, 0);
            
            % Complement the new MSB of dec_low, dec_up and dec_tag
            dec_low = bitxor(dec_low, 2^(N-1) );
            dec_up  = bitxor(dec_up,  2^(N-1) );
            dec_tag = bitxor(dec_tag, 2^(N-1) );

        end
    end % end while
end % end while length(dseq)

% Set the same output orientation as code
if (row_cd > 1)
    dseq = dseq.';
end
%-------------------------------------------------------------
function [ptr] = pick(cum_counts, value);
% This internal function is used to find where value is positioned

% Check for this case and quickly exit
if value == cum_counts(end)
    ptr = length(cum_counts)-1;
    return
end

c = find(cum_counts <= value);
ptr = c(end);

% EOF


function b = de2bi(varargin) 
%DE2BI Convert decimal numbers to binary numbers. 
%   B = DE2BI(D) converts a nonnegative integer decimal vector D to a binary 
%   matrix B. Each row of the binary matrix B corresponds to one element of D. 
%   The default orientation of the of the binary output is Right-MSB; the first 
%   element in B represents the lowest bit. 
% 
%   In addition to the vector input, three optional parameters can be given: 
% 
%   B = DE2BI(...,N) uses N to define how many digits (columns) are output. 
% 
%   B = DE2BI(...,N,P) uses P to define which base to convert the decimal 
%   elements to. 
% 
%   B = DE2BI(...,FLAG) uses FLAG to determine the output orientation.  FLAG 
%   has two possible values, 'right-msb' and 'left-msb'.  Giving a 'right-msb' 
%   FLAG does not change the function's default behavior.  Giving a 'left-msb' 
%   FLAG flips the output orientation to display the MSB to the left. 
% 
%   Examples: 
%   » D = [12; 5]; 
% 
%   » B = de2bi(D)                  » B = de2bi(D,5) 
%   B =                             B = 
%        0     0     1     1             0     0     1     1     0 
%        1     0     1     0             1     0     1     0     0 
% 
%   » T = de2bi(D,[],3)             » B = de2bi(D,5,'left-msb') 
%   T =                             B = 
%        0     1     1                   0     1     1     0     0 
%        2     1     0                   0     0     1     0     1 
% 
%   See also BI2DE. 
 
%   Copyright 1996-2001 The MathWorks, Inc. 
%   $Revision: 1.16 $  $Date: 2001/04/23 15:32:11 $ 
 
% Typical error checking. 
error(nargchk(1,4,nargin)); 
 
% --- Placeholder for the signature string. 
sigStr = ''; 
flag = ''; 
p = []; 
n = []; 
 
% --- Identify string and numeric arguments 
for i=1:nargin 
   if(i>1) 
      sigStr(size(sigStr,2)+1) = '/'; 
   end; 
   % --- Assign the string and numeric flags 
   if(ischar(varargin{i})) 
      sigStr(size(sigStr,2)+1) = 's'; 
   elseif(isnumeric(varargin{i})) 
      sigStr(size(sigStr,2)+1) = 'n'; 
   else 
      error('Only string and numeric arguments are accepted.'); 
   end; 
end; 
 
% --- Identify parameter signitures and assign values to variables 
switch sigStr 
   % --- de2bi(d) 
   case 'n' 
      d		= varargin{1}; 
 
	% --- de2bi(d, n) 
	case 'n/n' 
      d		= varargin{1}; 
      n		= varargin{2}; 
 
	% --- de2bi(d, flag) 
	case 'n/s' 
      d		= varargin{1}; 
      flag	= varargin{2}; 
 
	% --- de2bi(d, n, flag) 
	case 'n/n/s' 
      d		= varargin{1}; 
      n		= varargin{2}; 
      flag	= varargin{3}; 
 
	% --- de2bi(d, flag, n) 
	case 'n/s/n' 
      d		= varargin{1}; 
      flag	= varargin{2}; 
      n		= varargin{3}; 
 
	% --- de2bi(d, n, p) 
	case 'n/n/n' 
      d		= varargin{1}; 
      n		= varargin{2}; 
      p  	= varargin{3}; 
 
	% --- de2bi(d, n, p, flag) 
	case 'n/n/n/s' 
      d		= varargin{1}; 
      n		= varargin{2}; 
      p  	= varargin{3}; 
      flag	= varargin{4}; 
 
	% --- de2bi(d, n, flag, p) 
	case 'n/n/s/n' 
      d		= varargin{1}; 
      n		= varargin{2}; 
      flag	= varargin{3}; 
      p  	= varargin{4}; 
 
	% --- de2bi(d, flag, n, p) 
	case 'n/s/n/n' 
      d		= varargin{1}; 
      flag	= varargin{2}; 
      n		= varargin{3}; 
      p  	= varargin{4}; 
 
   % --- If the parameter list does not match one of these signatures. 
   otherwise 
      error('Syntax error.'); 
end; 
 
if isempty(d) 
   error('Required parameter empty.'); 
end 
 
d = d(:); 
len_d = length(d); 
 
if max(max(d < 0)) | max(max(~isfinite(d))) | (~isreal(d)) | (max(max(floor(d) ~= d))) 
   error('Input must contain only finite real positive integers.'); 
end 
 
% Assign the base to convert to. 
if isempty(p) 
    p = 2; 
elseif max(size(p) ~= 1) 
   error('Destination base must be scalar.'); 
elseif (~isfinite(p)) | (~isreal(p)) | (floor(p) ~= p) 
   error('Destination base must be a finite real integer.'); 
elseif p < 2 
   error('Cannot convert to a base of less than two.'); 
end; 
 
% Determine minimum length required. 
tmp = max(d); 
if tmp ~= 0 				% Want base-p log of tmp. 
   ntmp = floor( log(tmp) / log(p) ) + 1; 
else 							% Since you can't take log(0). 
   ntmp = 1; 
end 
 
% This takes care of any round off error that occurs for really big inputs. 
if ~( (p^ntmp) > tmp ) 
   ntmp = ntmp + 1; 
end 
 
% Assign number of columns in output matrix. 
if isempty(n) 
   n = ntmp; 
elseif max(size(n) ~= 1) 
   error('Specified number of columns must be scalar.'); 
elseif (~isfinite(n)) | (~isreal(n)) | (floor(n) ~= n) 
   error('Specified number of columns must be a finite real integer.'); 
elseif n < ntmp 
   error('Specified number of columns in output matrix is too small.'); 
end 
 
% Check if the string flag is valid. 
if isempty(flag) 
   flag = 'right-msb'; 
elseif ~(strcmp(flag, 'right-msb') | strcmp(flag, 'left-msb')) 
   error('Invalid string flag.'); 
end 
 
% Initial value. 
b = zeros(len_d, n); 
 
% Perform conversion. 
for i = 1 : len_d                   % Cycle through each element of the input vector/matrix. 
    j = 1; 
    tmp = d(i); 
    while (j <= n) & (tmp > 0)      % Cycle through each digit. 
        b(i, j) = rem(tmp, p);      % Determine current digit. 
        tmp = floor(tmp/p); 
        j = j + 1; 
    end; 
end; 
 
% If a flag is specified to flip the output such that the MSB is to the left. 
if strcmp(flag, 'left-msb') 
 
   b2 = b; 
   b = b2(:,n:-1:1); 
 
end 
 
% [EOF] de2bi.m 

function d = bi2de(b, varargin) 
%BI2DE Convert binary vectors to decimal numbers. 
%   D = BI2DE(B) converts a binary vector B to a decimal value D. When B is a 
%   matrix, the conversion is performed row-wise and the output D is a column 
%   vector of decimal values. The default orientation of the binary input 
%   is Right-MSB; the first element in B represents the least significant bit. 
% 
%   In addition to the input matrix, two optional parameters can be given: 
% 
%   D = BI2DE(...,P) converts a base P vector to a decimal value. 
% 
%   D = BI2DE(...,FLAG) uses FLAG to determine the input orientation.  FLAG has 
%   two possible values, 'right-msb' and 'left-msb'.  Giving a 'right-msb' FLAG 
%   does not change the function's default behavior.  Giving a 'left-msb' 
%   FLAG flips the input orientation such that the MSB is on the left. 
% 
%   Examples: 
%   » B = [0 0 1 1; 1 0 1 0]; 
%   » T = [0 1 1; 2 1 0]; 
% 
%   » D = bi2de(B)      » D = bi2de(B,'left-msb')      » D = bi2de(T,3) 
%   D =                 D =                            D = 
%       12                   3                             12 
%        5                  10                              5 
% 
%   See also DE2BI. 
 
%   Copyright 1996-2004 The MathWorks, Inc. 
%   $Revision: 1.15.4.1 $  $Date: 2004/08/10 01:32:22 $ 
 
 
% Typical error checking. 
error(nargchk(1,3,nargin)); 
 
% --- Placeholder for the signature string. 
sigStr = ''; 
flag = ''; 
p = []; 
 
% Check the type of the input B 
if ~(isnumeric(b) || islogical(b)) 
    error('The binary input must be numeric or logical.'); 
end 
b = double(b);  % To allow logicals to work 
 
% --- Identify string and numeric arguments 
for i=1:length(varargin) 
   if(i>1) 
      sigStr(size(sigStr,2)+1) = '/'; 
   end 
   % --- Assign the string and numeric flags 
   if(ischar(varargin{i})) 
      sigStr(size(sigStr,2)+1) = 's'; 
   elseif(isnumeric(varargin{i})) 
      sigStr(size(sigStr,2)+1) = 'n'; 
   else 
      error('Optional parameters must be string or numeric.'); 
   end 
end 
 
% --- Identify parameter signitures and assign values to variables 
switch sigStr 
     
    % --- bi2de(d) 
    case '' 
         
	% --- bi2de(d, p) 
	case 'n' 
      p		= varargin{1}; 
 
	% --- bi2de(d, flag) 
	case 's' 
      flag	= varargin{1}; 
 
	% --- bi2de(d, p, flag) 
	case 'n/s' 
      p		= varargin{1}; 
      flag	= varargin{2}; 
 
	% --- bi2de(d, flag, p) 
	case 's/n' 
      flag	= varargin{1}; 
      p		= varargin{2}; 
 
   % --- If the parameter list does not match one of these signatures. 
   otherwise 
      error('Syntax error.'); 
end 
 
if isempty(b) 
   error('Required parameter empty.'); 
end 
 
if max(max(b < 0)) || max(max(~isfinite(b))) || (~isreal(b)) || ... 
     (max(max(floor(b) ~= b))) 
    error('Input must contain only finite real positive integers.'); 
end 
 
% Set up the base to convert from. 
if isempty(p) 
    p = 2; 
elseif max(size(p)) > 1 
   error('Source base must be a scalar.'); 
elseif (floor(p) ~= p) || (~isfinite(p)) || (~isreal(p)) 
   error('Source base must be a finite real integer.'); 
elseif p < 2 
   error('Source base must be greater than or equal to two.'); 
end 
 
if max(max(b)) > (p-1) 
   error('The elements of the matrix are larger than the base can represent.'); 
end 
 
n = size(b,2); 
 
% If a flag is specified to flip the input such that the MSB is to the left. 
if isempty(flag) 
   flag = 'right-msb'; 
elseif ~(strcmp(flag, 'right-msb') || strcmp(flag, 'left-msb')) 
   error('Invalid string flag.'); 
end 
 
if strcmp(flag, 'left-msb') 
 
   b2 = b; 
   b = b2(:,n:-1:1); 
 
end 
 
%%% The conversion 
max_length = 1024; 
pow2vector = p.^(0:1:(size(b,2)-1)); 
size_B = min(max_length,size(b,2)); 
d = b(:,1:size_B)*pow2vector(:,1:size_B).'; 
 
% handle the infs... 
idx = find(max(b(:,max_length+1:size(b,2)).') == 1); 
d(idx) = inf; 
 
% [EOF] 
