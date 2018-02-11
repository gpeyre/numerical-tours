function perform_toolbox_installation(varargin)

// perform_toolbox_installation - add toolboxes to the path
//
//   perform_toolbox_installation('signal', 'graph', 'general');
//
//   Copyright (c) 2010 Gabriel Peyre

for i=1:argn(2)
    getd(['toolbox_' varargin(i) '/']);
end