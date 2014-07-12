k = 100;
Tlist = linspace(0, 1, k);
err = [];
for i=1:k
    x1 = x .* (abs(x)>Tlist(i));
    err(i) = norm( x0-x1, 'fro' );
end
clf;
plot(Tlist, err); axis('tight'); 
set_label('T', 'Error');
