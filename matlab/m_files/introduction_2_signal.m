%% Introduction to Signal Processing
% This numerical tour explores some basic signal processing tasks.


perform_toolbox_installation('signal', 'general');

%% Loading and Displaying Signals
% Signals are 1D vectors, usually stored as |(n,1)| arrays, where |n| is the number of samples.

n = 512;

%%
% Load a signal. (function load_signal.m should be in the toolbox of each course)

f = load_signal('Piece-Regular', n); % signal of size n

%% 
% One can force to be a column vector (just to be sure).

f = f(:); 

%%
% One can rescale to [0,1] the entries of the signal.

f = rescale(f);

%% 
% Display the signal.

clf;
plot(1:n, f);
axis('tight');
title('My title'); % title
set_label('variable x', 'variable y'); % axis

%% 
% You can display several figures using |subplot|

% divide the screen in 2x2 and select 1st quadrant
subplot(2, 2, 1);
plot(f); axis('tight');
% select the last quadrant
subplot(2, 2, 4);
plot(f.^2); axis('tight');

%%
% You can display several signals on the same figure

clf;
plot(1:n, [f f.^2]');
legend('signal', 'signal^2');