function r = myisfield(options, f)

// myisfield - emulation of isfield
//
//	r = myisfield(options, f);
//
//	test if options.f exists.
//
//	Copyright (c) 2008 Gabriel Peyre

r = getfield(1,options);
r = or(r(3:$)==f);

endfunction