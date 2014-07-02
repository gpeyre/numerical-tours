function [fx,fy,fz] = grad(M, options)

% grad - gradient, forward differences
%
%   [gx,gy] = grad(M, options);
% or
%   g = grad(M, options);
%
%   options.bound = 'per' or 'sym'
%   options.order = 1 (backward differences)
%                 = 2 (centered differences)
%
%   Works also for 3D array.
%   Assme that the function is evenly sampled with sampling step 1.
%
%   See also: div.
%
%   Copyright (c) Gabriel Peyre


options.null = 0;
bound = getoptions(options, 'bound', 'sym');
order = getoptions(options, 'order', 1);


% retrieve number of dimensions
nbdims = 2;
if size(M,1)==1 || size(M,2)==1
    nbdims = 1;
end
if size(M,1)>1 && size(M,2)>1 && size(M,3)>1
    nbdims = 3;
end


if strcmp(bound, 'sym')    
    if order==1
        fx = M([2:end end],:,:)-M;
    else
        fx = ( M([2:end end],:,:)-M([1 1:end-1],:,:) )/2;
        % boundary
        fx(1,:,:) = M(2,:,:)-M(1,:,:);
        fx(end,:,:) = M(end,:,:)-M(end-1,:,:);
    end
    if nbdims>=2
        if order==1
            fy = M(:,[2:end end],:)-M;
        else
            fy = ( M(:,[2:end end],:)-M(:,[1 1:end-1],:) )/2;
            % boundary
            fy(:,1,:) = M(:,2,:)-M(:,1,:);
            fy(:,end,:) = M(:,end,:)-M(:,end-1,:);
        end
    end
    if nbdims>=3
        if order==1
            fz = M(:,:,[2:end end])-M;
        else
            fz = ( M(:,:,[2:end end])-M(:,:,[1 1:end-1]) )/2;
            % boundary
            fz(:,:,1) = M(:,:,2)-M(:,:,1);
            fz(:,:,end) = M(:,:,end)-M(:,:,end-1);
        end
    end
else
    if order==1
        fx = M([2:end 1],:,:)-M;
    else
        fx = ( M([2:end 1],:,:)-M([end 1:end-1],:,:) )/2;
    end
    if nbdims>=2
        if order==1
            fy = M(:,[2:end 1],:)-M;
        else
            fy = ( M(:,[2:end 1],:)-M(:,[end 1:end-1],:) )/2;
        end
    end
    if nbdims>=3
        if order==1
            fz = M(:,:,[2:end 1])-M;
        else
            fz = ( M(:,:,[2:end 1])-M(:,:,[end 1:end-1]) )/2;
        end
    end
end

if nargout==1
    if nbdims==2
        fx = cat(3,fx,fy);
    elseif nbdims==3
        fx = cat(4,fx,fy,fz);
    end
end