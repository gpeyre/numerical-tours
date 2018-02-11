function r = isfield(options, f)

// isfield - emulation of isfield
//
//	r = isfield(options, f);
//
//	test if options.f exists.
//
//	Copyright (c) 2008 Gabriel Peyre

r = getfield(1,options);
r = or(r(3:$)==f);

endfunction