criter = [];
dlist = 1:N/20;
for i=1:length(dlist)
    s = twosparse(dlist(i));
    I = supp(s);
    criter(i,:) = [F(Phi,s) erc(Phi,I) werc(Phi,I)];        
end
criter(criter<0) = Inf;
clf; hold on;
h = plot(dlist, criter); set(h, 'LineWidth', 2);
h = plot(dlist, dlist*0+1, 'k--'); set(h, 'LineWidth', 2);
xlabel('d');  legend('F', 'ERC', 'w-ERC');
axis([1 max(dlist) min(criter(:)) 1.5]);
axis tight;
