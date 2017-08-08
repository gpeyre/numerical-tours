function extend_stack_size(mu, old_version)

// extend_stack_size - extend memory
//
//  Call
//      extend_stack_size();
//  to set stacksize to maximum.
//
//  If you have an old Scilab version and the call 
//  to extend_stack_size gives an error, then use
//
//  extend_stack_size(mult,1);
//
//  with mult a multiplying factor (e.g. 4).
//
//  Copyright (c) 2008 Gabriel Peyre


if argn(2)==1
    u = getversion();
    u = str2code(u);
    old_version = u(8)<5;
end


if old_version==0
    stacksize('max');
    return;
end

global stacksize_global;
global stacksize_max_global;
if isempty(stacksize_global)
    stacksize_global = 1;
end
if isempty(stacksize_max_global)
    stacksize_max_global = 20;
end

if argn(2)<1
    mu = 4;
end

if stacksize_global>stacksize_max_global
    warning(strcat(['Stacksize multiplier extend ' num2str(stacksize_max_global) '.']));
    warning(strcat(['Solution: do not call extend_stack_size too many time or extend global variable stacksize_max_global.']));
    return;
end

stacksize_global = stacksize_global*mu;

a = stacksize();
stacksize(mu*a(1));

endfunction