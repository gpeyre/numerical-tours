cR = sort(abs(fF(:)));
if cR(n^2)>cR(1)
    cR = reverse(cR); % be sure it is in reverse order
end
lw = 2;
clf;
h = plot(log10(cR)); 
if using_matlab()  
    set(h, 'LineWidth', lw);
end
axis('tight');
