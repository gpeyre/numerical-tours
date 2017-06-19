module nt_general

using PyPlot

export rescale, reverse, clamp

#####################################
function rescale(x,a=0,b=1)

	# rescale - rescale data in [a,b]
	#
	#   y = rescale(x,a,b);
	#
	#   Copyright (c) 2004 Gabriel Peyre

	m = minimum(x[:]);
	M = maximum(x[:]);

	if M-m<1e-10
	    y = x;
	else
	    y = (b-a) * (x-m)/(M-m) + a;
	end
	return y;

end ## end of rescale


#####################################
function reverse(x)
	# flip a vector
	return x[length(x):-1:1];
end  # function x = reverse(x)

#####################################
function clamp(x,a=0,b=1)

	# clamp - clamp a value
	#
	#   y = clamp(x,a,b);
	#
	# Default is [a,b]=[0,1].
	#
	#   Copyright (c) 2004 Gabriel Peyré

	y = max(x,a);
	y = min(y,b);
	return y;

end # clamp

end ## end of module