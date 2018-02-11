function set_rand_seeds(a,b)

// set_rand_seeds - initialize rand and randn

if argn(1)<1
    a = 123456;
end
if argn(2)<1
    b = 123456;
end
rand('normal');
randn('state', a);
rand('uniform');
randn('state', a);

endfunction