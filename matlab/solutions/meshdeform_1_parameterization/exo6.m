p = length(B);
t = linspace(0,2*pi(),p+1)'; t(p) = [];
x0 = cos(t); y0 = sin(t);
% 
Rx = zeros(n,1); Rx(B) = x0;
Ry = zeros(n,1); Ry(B) = y0;
if using_matlab()
    x = L1 \ Rx; 
    y = L1 \ Ry;
else
    options.maxit = 300;
    x = perform_cg(L1,Rx,options);
    y = perform_cg(L1,Ry,options);
end
Y = [x';y'];
clf;
options.lighting = 0;
plot_mesh([Y;zeros(1,n)],F);
shading('faceted');
