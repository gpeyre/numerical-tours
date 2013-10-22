function [selx,sely] = compute_quadsel(j,q, options)

% compute_quadsel - compute the indices for selecting subband of a wavelet transform.
%
%   [selx,sely] = compute_quadrant_selection(j,q, options);
%
%   q=1 is for phi/psi (on X / on Y)
%   q=2 is for psi/phi (on X / on Y)
%   q=3 is for psi/psi (on X / on Y)
%
%   % Suppose that MW contains a wavelet transform 
%   j = 4; q = 2;
%   [selx,sely] = compute_quadrant_selection(j,q);
%   MWs = MW(selx,sely);    % MWs contains the quadrant.
%
%   Depending on the transform you are using, you should 
%   set options.transform to 
%       'isotropic2d' (default): 2D standard wavelet basis.
%       'quincunx': for quincunx transform (perform_quincunx_wavelet_transform).
%       
%
%   Copyright (c) 2004 Gabriel Peyre

options.null = 0;
if isfield(options, 'transform')
    transform = options.transform;
else
    transform = 'isotropic2d';
end

switch lower(transform)
    case 'isotropic2d'
        if q<0 || q>3
            error('Quadrant number is in {0,1,2,3}.');
        end
        if q==0
            selx = 1:2^j;
            sely = 1:2^j;
        elseif q==1 % 1st quadrant
            selx = 1:2^j;
            sely = (2^j+1):2^(j+1);
        elseif q==2
            selx = (2^j+1):2^(j+1);
            sely = 1:2^j;
        elseif q==3
            selx = (2^j+1):2^(j+1);
            sely = (2^j+1):2^(j+1);
        else
            selx = []; sely = [];
        end
    case 'quincunx'
        if q==0
            selx = 1:2^j;
            sely = 1:2^j;            
        elseif q==1
            selx = 1:2^(j+1);
            sely = 2^j+1:2^(j+1);
        elseif q==2
            selx = 2^j+1:2^(j+1);
            sely = 1:2^j;
        else
            selx = []; sely = [];
        end
    otherwise 
        error('Unknown transform.');
end