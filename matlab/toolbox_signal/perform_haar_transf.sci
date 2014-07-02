function f = perform_haar_transf(f, Jmin, direction, options)

// perform_haar_transf - peform fast Haar transform
//
//   y = perform_haar_transf(x, Jmin, direction);
//
//   Implement a Haar wavelets.
//   Works in any dimension.
//
//   Copyright (c) 2008 Gabriel Peyre

n = size(f,1); 
Jmax = log2(n)-1; 

if direction==1
    ////// FORWARD //////
    for j=Jmax:-1:Jmin
        sel = 1:2^(j+1);
        a = subselect(f,sel);
        for d=1:nb_dims(f)
            Coarse = ( subselectdim(a,1:2:size(a,d),d) + subselectdim(a,2:2:size(a,d),d) )/sqrt(2);
            Detail = ( subselectdim(a,1:2:size(a,d),d) - subselectdim(a,2:2:size(a,d),d) )/sqrt(2);
            a = cat3(d, Coarse, Detail );
        end
        f = subassign(f,sel,a);
    end
else
    ////// BACKWARD //////
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

//////////////////////////////////////////////////////////////////////////////////
function f = subassigndim(f,sel,g,d)


if d==1
        f(sel,:,:,:,:,:,:,:) = g;
elseif d==2
        f(:,sel,:,:,:,:,:,:) = g;
elseif d==3
        f(:,:,sel,:,:,:,:,:) = g;
elseif d==4
        f(:,:,:,sel,:,:,:,:) = g;
elseif d==5
        f(:,:,:,:,sel,:,:,:) = g;
elseif d==6
        f(:,:,:,:,:,sel,:,:) = g;
elseif d==7
        f(:,:,:,:,:,:,sel,:) = g;
elseif d==8
        f(:,:,:,:,:,:,:,sel) = g;
else
        error('Not implemented');
end

endfunction