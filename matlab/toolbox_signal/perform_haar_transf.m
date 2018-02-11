function f = perform_haar_transf(f, Jmin, dir, options)

% perform_haar_transf - peform fast Haar transform
%
%   y = perform_haar_transf(x, Jmin, dir);
%
%   Implement a Haar wavelets.
%   Works in any dimension.
%
%   Copyright (c) 2008 Gabriel Peyre

n = size(f,1); 
Jmax = log2(n)-1; 

if dir==1
    %%% FORWARD %%%
    for j=Jmax:-1:Jmin
        sel = 1:2^(j+1);
        a = subselect(f,sel);
        for d=1:nb_dims(f)
            Coarse = ( subselectdim(a,1:2:size(a,d),d) + subselectdim(a,2:2:size(a,d),d) )/sqrt(2);
            Detail = ( subselectdim(a,1:2:size(a,d),d) - subselectdim(a,2:2:size(a,d),d) )/sqrt(2);
            a = cat(d, Coarse, Detail );
        end
        f = subassign(f,sel,a);
    end
else
    %%% BACKWARD %%%
    for j=Jmin:Jmax
        sel = 1:2^(j+1);
        a = subselect(f,sel);
        for d=1:nb_dims(f)
            Detail = subselectdim(a,2^j+1:2^(j+1),d);
            Coarse = subselectdim(a,1:2^j,d);
            a = subassigndim(a, 1:2:2^(j+1), ( Coarse + Detail )/sqrt(2),d );
            a = subassigndim(a, 2:2:2^(j+1), ( Coarse - Detail )/sqrt(2),d );
        end
        f = subassign(f,sel,a);
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subselect(f,sel)
switch nb_dims(f)
    case 1
        f = f(sel);
    case 2
        f = f(sel,sel);
    case 3
        f = f(sel,sel,sel);
    case 4
        f = f(sel,sel,sel,sel);
    case 5
        f = f(sel,sel,sel,sel,sel);
    case 6
        f = f(sel,sel,sel,sel,sel,sel);
    case 7
        f = f(sel,sel,sel,sel,sel,sel,sel);
    case 8
        f = f(sel,sel,sel,sel,sel,sel,sel,sel);
    otherwise
        error('Not implemented');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subselectdim(f,sel,d)
switch d
    case 1
        f = f(sel,:,:,:,:,:,:,:);
    case 2
        f = f(:,sel,:,:,:,:,:,:);
    case 3
        f = f(:,:,sel,:,:,:,:,:);
    case 4
        f = f(:,:,:,sel,:,:,:,:);
    case 5
        f = f(:,:,:,:,sel,:,:,:);
    case 6
        f = f(:,:,:,:,:,sel,:,:);
    case 7
        f = f(:,:,:,:,:,:,sel,:);
    case 8
        f = f(:,:,:,:,:,:,:,sel);
    otherwise
        error('Not implemented');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subassign(f,sel,g)
switch nb_dims(f)
    case 1
        f(sel) = g;
    case 2
        f(sel,sel) = g;
    case 3
        f(sel,sel,sel) = g;
    case 4
        f(sel,sel,sel,sel) = g;
    case 5
        f(sel,sel,sel,sel,sel) = g;
    case 6
        f(sel,sel,sel,sel,sel,sel) = g;
    case 7
        f(sel,sel,sel,sel,sel,sel,sel) = g;
    case 8
        f(sel,sel,sel,sel,sel,sel,sel,sel) = g;
    otherwise
        error('Not implemented');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = subassigndim(f,sel,g,d)
switch d
    case 1
        f(sel,:,:,:,:,:,:,:) = g;
    case 2
        f(:,sel,:,:,:,:,:,:) = g;
    case 3
        f(:,:,sel,:,:,:,:,:) = g;
    case 4
        f(:,:,:,sel,:,:,:,:) = g;
    case 5
        f(:,:,:,:,sel,:,:,:) = g;
    case 6
        f(:,:,:,:,:,sel,:,:) = g;
    case 7
        f(:,:,:,:,:,:,sel,:) = g;
    case 8
        f(:,:,:,:,:,:,:,sel) = g;
    otherwise
        error('Not implemented');
end