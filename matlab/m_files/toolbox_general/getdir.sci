function getdir(dname)

// getdir - replacement for getd
// 
//  getdir(dirname);
//
//  Copyright (c) 2008 Gabriel Peyre

a = dir(strcat([dname '*.sci']));

for i=1:size(a.name, 1)
    filename = a.name(i);
    exec(filename);
 //   exec("toolbox_general/atan2.sci");
end

endfunction
