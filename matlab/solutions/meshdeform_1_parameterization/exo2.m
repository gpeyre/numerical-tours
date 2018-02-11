R = zeros(2, n);
R(:,B) = Z;
if using_matlab()
    Y = (L1 \ R')'; 
else
    options.maxit = 300;
    Y(1,:) = perform_cg(L1,R(1,:)',options)';
    Y(2,:) = perform_cg(L1,R(2,:)',options)';
end
clf;
plot_mesh([Y;zeros(1,n)],F);
shading faceted; axis tight;
