function f = perform_wavortho_transf(f,Jmin,dir,options)

% perform_wavortho_transf - compute orthogonal wavelet transform
%
%   fw = perform_wavortho_transf(f,Jmin,dir,options);
%
%   You can give the filter in options.h.
%
%   Works in arbitrary dimension.
%
%   Copyright (c) 2009 Gabriel Peyre

options.null = 0;
h = getoptions(options,'h', compute_wavelet_filter('Daubechies',4) );
g = [0 h(length(h):-1:2)] .* (-1).^(1:length(h));

n = size(f,1); 
Jmax = log2(n)-1; 

if dir==1
    %%% FORWARD %%%
    for j=Jmax:-1:Jmin
        sel = 1:2^(j+1);
        a = subselect(f,sel);
        for d=1:nb_dims(f)
            a = cat(d, subsampling(cconv(a,h,d),d), subsampling(cconv(a,g,d),d) );
        end
        f = subassign(f,sel,a);
    end
else
    %%% FORWARD %%%
    for j=Jmin:Jmax
        sel = 1:2^(j+1);
        a = subselect(f,sel);
        for d=1:nb_dims(f)
            w = subselectdim(a,2^j+1:2^(j+1),d);
            a = subselectdim(a,1:2^j,d);
            a = cconv(upsampling(a,d),reverse(h),d) + cconv(upsampling(w,d),reverse(g),d);
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