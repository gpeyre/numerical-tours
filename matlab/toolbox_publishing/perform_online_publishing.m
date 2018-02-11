function perform_online_publishing(name,location)

% perform_online_publishing - copy a numerical tour to its final destination
%
%   perform_online_publishing(name,location);
%
%   Copyright (c) 2010 Gabriel Peyre

if nargin<2
    location = '/Volumes/Ceremade/';
end

% copy the files
src = ['../html/' name '/'];
tgt = [location 'numerical-tour/tours/' name '/'];
system(['mkdir ' tgt]);
system(['cp ' src '* ' tgt]);

% copy the index files
perform_index_generation();
src = ['../html/index_tours.php'];
tgt = [location 'numerical-tour/tours/'];
system(['cp ' src ' ' tgt]);
