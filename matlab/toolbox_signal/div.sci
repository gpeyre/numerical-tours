function fd = div(Px,Py, options)

// div - divergence operator
//
//	fd = div(Px,Py, options);
//	fd = div(P, options);
//
//   options.bound = 'per' or 'sym'
//   options.order = 1 (backward differences)
//                 = 2 (centered differences)
//
//   Note that the div and grad operator are adjoint
//   of each other such that 
//       <grad(f),g>=<f,div(g)>
//
//   See also: grad.
//
//	Copyright (c) 2007 Gabriel Peyre

// retrieve number of dimensions
nbdims = 2;
if size(Px,1)==1 | size(Px,2)==1
    nbdims = 1;
end
if size(Px,3)>1
	if argn(2)>1
		options = Py;
        clear Py;
    end
        Py = Px(:,:,2);
        Px = Px(:,:,1);
//    if nd_dims(Px)==4 size(Px,4)<=1
//        Py = Px(:,:,2);
//        Px = Px(:,:,1);
//    else
//        Pz = Px(:,:,:,3);
//        Py = Px(:,:,:,2);
//        Px = Px(:,:,:,1); 
//        nbdims = 3;
//    end
end

options.null = 0;
bound = getoptions(options, 'bound', 'sym');
order = getoptions(options, 'order', 1);

if strcmp(bound, 'sym')
    if order==1
        fx = Px-Px([1 1:n-1],:,:);         
        fx(1,:,:)   = Px(1,:,:);        // boundary
        fx(n,:,:) = -Px(n-1,:,:);        
        if nbdims>=2
            fy = Py-Py(:,[1 1:n-1],:);
            fy(:,1,:)   = Py(:,1,:);    // boundary
            fy(:,n,:) = -Py(:,n-1,:);
        end        
        if nbdims>=3
            fz = Pz-Pz(:,:,[1 1:n-1]);  
            fz(:,:,1)   = Pz(:,:,1);    // boundary
            fz(:,:,n) = -Pz(:,:,n-1);
        end
    else
        fx = (Px([2:n n],:,:)-Px([1 1:n-1],:,:))/2;
        fx(1,:,:)   = +Px(2,:,:)/2+Px(1,:,:);   // boundary
        fx(2,:,:)   = +Px(3,:,:)/2-Px(1,:,:);
        fx(n,:,:) = -Px(n,:,:)-Px(n-1,:,:)/2;
        fx(n-1,:,:) = Px(n,:,:)-Px(n-2,:,:)/2;
        if nbdims>=2
            fy = (Py(:,[2:n n],:)-Py(:,[1 1:n-1],:))/2;
            fy(:,1,:)   = +Py(:,2,:)/2+Py(:,1,:);
            fy(:,2,:)   = +Py(:,3,:)/2-Py(:,1,:);
            fy(:,n,:) = -Py(:,n,:)-Py(:,n-1,:)/2;
            fy(:,n-1,:) = Py(:,n,:)-Py(:,n-2,:)/2;
        end
        if nbdims>=3
            fz = (Pz(:,:,[2:n n])-Pz(:,:,[1 1:n-1]))/2;            
            fz(:,:,1)   = +Pz(:,:,2)/2+Pz(:,:,1);   // boundary
            fz(:,:,2)   = +Pz(:,:,3)/2-Pz(:,:,1);
            fz(:,:,n) = -Pz(:,:,n)-Pz(:,:,n-1)/2;
            fz(:,:,n-1) = Pz(:,:,n)-Pz(:,:,n-2)/2;
        end
    end 
else
    if order==1
        fx = Px-Px([n 1:n-1],:,:);
        if nbdims>=2
            fy = Py-Py(:,[n 1:n-1],:);
        end
        if nbdims>=3
            fz = Pz-Pz(:,:,[n 1:n-1]);
        end
    else
        fx = (Px([2:n 1],:,:)-Px([n 1:n-1],:,:))/2;
        if nbdims>=2
            fy = (Py(:,[2:n 1],:)-Py(:,[n 1:n-1],:))/2;
        end
        if nbdims>=3
            fz = (Pz(:,:,[2:n 1])-Pz(:,:,[n 1:n-1]))/2;
        end
    end
end

// gather result
if nbdims==3
    fd = fx+fy+fz;
elseif nbdims==2
    fd = fx+fy;
else
    fd = fx;
end

endfunction