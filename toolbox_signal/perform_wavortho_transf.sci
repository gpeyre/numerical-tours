function f = perform_wavortho_transf(f,Jmin,direction,options)

// perform_wavortho_transf - compute orthogonal wavelet transform
//
//   fw = perform_wavortho_transf(f,Jmin,direction,options);
//
//   You can give the filter in options.h.
//
//   Works in arbitrary dimension.
//
//   Copyright (c) 2009 Gabriel Peyre

options.null = 0;
h = getoptions(options,'h', compute_wavelet_filter('Daubechies',4) );
g = [0 h(length(h):-1:2)] .* (-1).^(1:length(h));

n = size(f,1); 
Jmax = log2(n)-1; 

if direction==1
    ////// FORWARD //////
    for j=Jmax:-1:Jmin
        sel = 1:2^(j+1);
        a = subselect(f,sel);
        for d=1:nb_dims(f)
            a = cat3(d, subsampling(cconv(a,h,d),d), subsampling(cconv(a,g,d),d) );
        end
        f = subassign(f,sel,a);
    end
else
    ////// FORWARD //////
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

endfunction
//////////////////////////////////////////////////////////////////////////////////
function f = subselect(f,sel)

d = nb_dims(f);

if d==1
        f = f(sel);
elseif d==2
        f = f(sel,sel);
elseif d==3
        f = f(sel,sel,sel);
elseif d==4
        f = f(sel,sel,sel,sel);
elseif d==5
        f = f(sel,sel,sel,sel,sel);
elseif d==6
        f = f(sel,sel,sel,sel,sel,sel);
elseif d==7
        f = f(sel,sel,sel,sel,sel,sel,sel);
elseif d==8
        f = f(sel,sel,sel,sel,sel,sel,sel,sel);
else
        error('Not implemented');
end

endfunction
//////////////////////////////////////////////////////////////////////////////////
function f = subselectdim(f,sel,d)

if d==1
        f = f(sel,:,:,:,:,:,:,:);
elseif d==2
        f = f(:,sel,:,:,:,:,:,:);
elseif d==3
        f = f(:,:,sel,:,:,:,:,:);
elseif d==4
        f = f(:,:,:,sel,:,:,:,:);
elseif d==5
        f = f(:,:,:,:,sel,:,:,:);
elseif d==6
        f = f(:,:,:,:,:,sel,:,:);
elseif d==7
        f = f(:,:,:,:,:,:,sel,:);
elseif d==8
        f = f(:,:,:,:,:,:,:,sel);
else
        error('Not implemented');
end

endfunction
//////////////////////////////////////////////////////////////////////////////////
function f = subassign(f,sel,g)

d = nb_dims(f);

if d==1
        f(sel) = g;
elseif d==2
        f(sel,sel) = g;
elseif d==3
        f(sel,sel,sel) = g;
elseif d==4
        f(sel,sel,sel,sel) = g;
elseif d==5
        f(sel,sel,sel,sel,sel) = g;
elseif d==6
        f(sel,sel,sel,sel,sel,sel) = g;
elseif d==7
        f(sel,sel,sel,sel,sel,sel,sel) = g;
elseif d==8
        f(sel,sel,sel,sel,sel,sel,sel,sel) = g;
else
        error('Not implemented');
end


endfunction