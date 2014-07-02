function compute_movie_file(A, filename, options)

% compute_movie_file - create an avi or gif file
%
%   compute_movie_file(A, filename, options);
%
%   A can be a 3D array or a cell array of 2D matrices.
%   Each A{k} or A(:,:,k) is one frame in the animation.
%
%   Output avi or gif file depending on the extension of the filename.
%   (you can also use options.format='avi' or 'gif')
%
%   Copyright (c) 2007 Gabriel Peyre


options.null = 0;
fps = getoptions(options, 'fps', 10);
format = getoptions(options, 'format', []);

if isempty(format)
    % find extension
    i = strfind(filename,'.');
    if isempty(i)
        error('Does not recognize extension, you should specify options.format.');
    end
    format = filename(i(end)+1:end);
end

if not(iscell(A)) && not(isstr(A))
    B = A; A = {};
    for i=1:size(B,3)
        A{end+1} = B(:,:,i);
    end
    clear B;
end


% close all;
% figure;

if iscell(A)
    nb = length(A);
else
    nb = length( dir([A '*.png']) );
end

if strcmp(format, 'avi')
    % 'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE' or 'None'
    compressor = 'Cinepak';
    compressor = 'none';
    qual = 50;
    aviobj = avifile(filename, 'compression', compressor, 'quality', qual, 'fps', fps );
    % aviobj = avifile(filename, 'fps', fps );
end

imap = getoptions(options, 'colormap', jet(256));
loopCount = getoptions(options, 'loopcount', 1);
verb = getoptions(options, 'verb', 0);

for j=1:nb
    if verb
        progressbar(j,nb);
    end
    if iscell(A)
        B = A{j};
    else
        % A is a filename
        B = load_image( [A '-' num2string_fixeddigit(j-1,3)] );
    end
    B = round(rescale( B, 1,256));
    switch format
        case 'avi'
            if size(B,3)==3
                frame = im2frame( rescale(B) );
            else
                frame = im2frame(B,gray(256));
            end
            % imageplot(A{j});
            % frame = getframe ( gca );
            aviobj = addframe( aviobj, frame );
        case 'gif'
            warning off;
            if j==1
                imwrite(B,imap,filename,'gif',...
                    'DelayTime',1/fps, 'LoopCount',loopCount,...
                    'WriteMode','overwrite')
            else
                imwrite(B,imap,filename,'gif','WriteMode','append',...
                    'DelayTime',1/fps);
            end
            warning on;
        otherwise
            error('Unknown format');
    end
end

if strcmp(format, 'avi')
    warning off;
    aviobj = close( aviobj );
    warning on;
end