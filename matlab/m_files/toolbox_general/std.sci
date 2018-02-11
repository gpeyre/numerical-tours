function y = std(x,a)

// std - for matlab compatibility
    

if argn(2)==1
    y = stdev(x);
else
    y = stdev(x,a);
end

endfunction