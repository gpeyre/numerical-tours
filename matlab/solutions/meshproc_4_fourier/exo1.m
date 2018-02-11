nblist = round( linspace(10,nb,4) );
clf;
for i=1:length(nblist)
    V = U(:,1:nblist(i));
    subplot(2,2,i);
    plot_mesh((vertex*V)*V',faces);
end
