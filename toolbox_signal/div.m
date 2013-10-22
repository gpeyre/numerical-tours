function fd = div(Px,Py, options)

% div - divergence operator
%
%	fd = div(Px,Py, options);
%	fd = div(P, options);
%
%   options.bound = 'per' or 'sym'
%   options.order = 1 (backward differences)
%                 = 2 (centered differences)
%
%   Note that the -div and grad operator are adjoint
%   of each other such that 
%       <grad(f),g>=<f,-div(g)>
%
%   See also: grad.
%
%	Copyright (c) 2007 Gabriel Peyre

% retrieve number of dimensions
nbdims = 2;
if size(Px,1)==1 || size(Px,2)==1
    nbdims = 1;
end
if size(Px,3)>1
	if nargin>1
		options = Py;
        clear Py;
    end
    if size(Px,4)<=1
        Py = Px(:,:,2);
        Px = Px(:,:,1);
    else
        Pz = Px(:,:,:,3);
        Py = Px(:,:,:,2);
        Px = Px(:,:,:,1); 
        nbdims = 3;
    end
end

options.null = 0;
bound = getoptions(options, 'bound', 'sym');
order = getoptions(options, 'order', 1);

if strcmp(bound, 'sym')
    if order==1
        fx = Px-Px([1 1:end-1],:,:);         
        fx(1,:,:)   = Px(1,:,:);        % boundary
        fx(end,:,:) = -Px(end-1,:,:);        
        if nbdims>=2
            fy = Py-Py(:,[1 1:end-1],:);
            fy(:,1,:)   = Py(:,1,:);    % boundary
            fy(:,end,:) = -Py(:,end-1,:);
        end        
        if nbdims>=3
            fz = Pz-Pz(:,:,[1 1:end-1]);  
            fz(:,:,1)   = Pz(:,:,1);    % boundary
            fz(:,:,end) = -Pz(:,:,end-1);
        end
    else
        fx = (Px([2:end end],:,:)-Px([1 1:end-1],:,:))/2;
        fx(1,:,:)   = +Px(2,:,:)/2+Px(1,:,:);   % boundary
        fx(2,:,:)   = +Px(3,:,:)/2-Px(1,:,:);
        fx(end,:,:) = -Px(end,:,:)-Px(end-1,:,:)/2;
        fx(end-1,:,:) = Px(end,:,:)-Px(end-2,:,:)/2;
        if nbdims>=2
            fy = (Py(:,[2:end end],:)-Py(:,[1 1:end-1],:))/2;
            fy(:,1,:)   = +Py(:,2,:)/2+Py(:,1,:);
            fy(:,2,:)   = +Py(:,3,:)/2-Py(:,1,:);
            fy(:,end,:) = -Py(:,end,:)-Py(:,end-1,:)/2;
            fy(:,end-1,:) = Py(:,end,:)-Py(:,end-2,:)/2;
        end
        if nbdims>=3
            fz = (Pz(:,:,[2:end end])-Pz(:,:,[1 1:end-1]))/2;            
            fz(:,:,1)   = +Pz(:,:,2)/2+Pz(:,:,1);   % boundary
            fz(:,:,2)   = +Pz(:,:,3)/2-Pz(:,:,1);
            fz(:,:,end) = -Pz(:,:,end)-Pz(:,:,end-1)/2;
            fz(:,:,end-1) = Pz(:,:,end)-Pz(:,:,end-2)/2;
        end
    end 
else
    if order==1
        fx = Px-Px([end 1:end-1],:,:);
        if nbdims>=2
            fy = Py-Py(:,[end 1:end-1],:);
        end
        if nbdims>=3
            fz = Pz-Pz(:,:,[end 1:end-1]);
        end
    else
        fx = (Px([2:end 1],:,:)-Px([end 1:end-1],:,:))/2;
        if nbdims>=2
            fy = (Py(:,[2:end 1],:)-Py(:,[end 1:end-1],:))/2;
        end
        if nbdims>=3
            fz = (Pz(:,:,[2:end 1])-Pz(:,:,[end 1:end-1]))/2;
        end
    end
end

% gather result
if nbdims==3
    fd = fx+fy+fz;
elseif nbdims==2
    fd = fx+fy;
else
    fd = fx;
end