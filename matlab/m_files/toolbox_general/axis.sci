function axis(a)

// axis - set axis bounds

if isstr(a)
	return;
end

g = get("current_axes");
g.data_bounds = a;

endfunction