t = [0 .02 .1 1];
LL0 = L0*L0;
clf;
for i=1:4
    LL1 = t(i)*L0 + (1-t(i))*LL0;
    LL1(I,:) = 0;
    LL1(I + (I-1)*n) = 1;
    % deform
    vertex1 = vertex0+( LL1 \ Delta0' )';
    % Display it.
    subplot(2,2,i);
    plot_mesh(vertex1,faces);
    view(-150,45);
    zoom(.8);
end
