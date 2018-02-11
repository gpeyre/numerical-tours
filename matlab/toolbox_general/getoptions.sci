function v = getoptions(options, name, v, mendatory)

// getoptions - retrieve options parameter
//
// v = getoptions(options, name, v, mendatory)
//
//  Copyright (c) 2008 Gabriel Peyre

if argn(2)<3
    error('Not enough arguments.');
end
if argn(2)<4
    mendatory = 0;
end

if isfield(options, name)
    v = eval( strcat(['options.' name ';']) );
elseif mendatory
    error( strcat(['You have to provide options.' name '.']) );
end

endfunction