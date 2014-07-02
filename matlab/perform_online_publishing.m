function perform_online_publishing(name,location)

% perform_online_publishing - copy a numerical tour to its final destination
%
%   perform_online_publishing(name,location);
%
%   mkdir /Volumes/Ceremade
%   echo Gwenn.26 | sshfs peyre@beta.ceremade.dauphine.fr:www/ /Volumes/Ceremade -o password_stdin,volname=Ceremade
%   unmount /Volumes/Ceremade
%   
%   Copyright (c) 2010 Gabriel Peyre

if nargin<2
    % mount location
    location = '/Volumes/Ceremade/';
end

% mount the file system
sshfs_cmd = '/usr/local/bin/sshfs';
system(['mkdir ' location]);
system(['echo U1AwAlFPqU | ' sshfs_cmd ' peyre@web.ceremade.dauphine.fr:www/ ' location ' -o password_stdin,volname=Ceremade']);

% copy the files
src = ['../html/' name '/'];
tgt = [location 'numerical-tour/tours/' name '/'];
system(['mkdir ' tgt]);
system(['cp ' src '* ' tgt]);

% copy the index files
perform_index_generation();
my_copy('../html/index_tours.php',  [location 'numerical-tour/tours/']);
my_copy('../html/index_news.php',   [location 'numerical-tour/tours/']);

% unmount the fs (not working)
% system(['umount ' location]);

end

%%%% AUX FUNCTIONS %%%
function my_copy(src, tgt)
system(['cp ' src ' ' tgt]);
end
