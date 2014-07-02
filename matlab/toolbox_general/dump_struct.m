function fid = dump_struct(s,fid,header)

% dump_struct - dump the content of a struct to a file
%
%   dump_struct(s,fid, header);
%
%   Copyright (c) 2008 Gabriel Peyre 

if nargin<3
    header = '';
end

if isstr(fid)
    fid = fopen(fid, 'a');
    if fid<=0
        error(['File ' fid ' does not exist.']);
    end
end

if not(isempty(header))
    fprintf(fid, '%s\n', header );
end

% central line (position of ':')
cl = 15;

if iscell(s)
    for i=1:length(s)
        fprintf(fid, '%s\n',  convert_string(s{i},0) );
    end    
elseif isstruct(s)
    s = orderfields(s);
    f = fieldnames(s);
    for i=1:length(f)
        v = getfield(s,f{i});
        fprintf(fid, '%s: %s\n',  convert_string(f{i},cl), convert_string(v,0) );
    end
else
    error('Unknown type.');
end

if nargout==0
    fclose(fid);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s = convert_string(v,cl)

if nargin<2
    cl = 15;
end

nmax = 10; % maximum of dumped string results

if isa(v, 'numeric')
    %    v = v(1);
    s = [];
    for i = 1:min(length(v(:)), nmax)
        if mod(v,1)==0
            s = [s int2str(v(i)) ' '];
        else
            s = [s num2str(v(i)) ' '];
        end
    end
    if length(v(:))>nmax
        s = [s '...'];
    end
elseif isstr(v)
    s = v;
else
    s = '-';
end

% try to align on central line
d = cl - length(s);
if d>0
    s = [repmat(' ',[1 d]) s];
end



