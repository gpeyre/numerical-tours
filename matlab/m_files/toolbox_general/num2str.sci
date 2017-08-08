function s = num2str(v,f)

// num2str - convert numerics to string
//
//	s = num2str(v);
//
//	(this function ensure matlab compatibility)
//
//	Copyright (c) 2008 Gabriel Peyre
	
if argn(2)<2
    f = 4;
end
s = mtlb_num2str(v);

endfunction